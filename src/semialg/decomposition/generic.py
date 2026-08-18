from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..cad.lifting.stack import CADCell
from ..domains import apply_assumptions, normalize_assumptions, normalize_domain
from ..formula import Formula, parse_formula, parse_formula_text, to_sympy
from ..generic import (
    ExceptionalCause,
    GenericSplit,
    generic_split_from_cells,
    input_boundary_causes,
    projection_causes,
)
from ..qe.complete import cells_to_formula
from ..simplify.result import simplify_qe_formula
from .cylindrical import (
    CADResult,
    CellSet,
    _cell_dim_with_levels,
    _formula_from_conditions,
    _normalize_formula,
    _normalize_variables,
    cad,
)

GenericOutput = Literal["formula", "cases", "cells", "function"]


@dataclass(frozen=True)
class GenericCase:
    """One parameter cell in a generic CAD analysis."""

    param_condition: sp.Expr
    solution_formula: sp.Expr
    param_cell: CADCell | None
    solution_cells: tuple[CADCell, ...]
    param_sample: Mapping[sp.Symbol, sp.Expr]
    dimension: int
    exceptional: bool = False

    @property
    def has_solution(self) -> bool:
        return self.solution_cells != () and self.solution_formula != sp.false

    def sample_point(self) -> Mapping[sp.Symbol, sp.Expr]:
        return self.param_sample

    def fiber_formula(self) -> sp.Expr:
        if not self.has_solution:
            return sp.false
        return simplify_qe_formula(
            self.solution_formula.subs(dict(self.param_sample)), implication_minimize=False
        )


@dataclass(frozen=True)
class GenericCADResult:
    """Result returned by :func:`generic_cad`."""

    generic_formula: sp.Expr
    exceptional_formula: sp.Expr
    cases: tuple[GenericCase, ...]
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    cad_result: CADResult
    output: GenericOutput
    status: str
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    generic_split: GenericSplit | None = None
    exceptional_causes: tuple[ExceptionalCause, ...] = ()

    @property
    def generic_cases(self) -> tuple[GenericCase, ...]:
        return tuple(case for case in self.cases if not case.exceptional and case.has_solution)

    @property
    def all_generic_cases(self) -> tuple[GenericCase, ...]:
        return tuple(case for case in self.cases if not case.exceptional)

    @property
    def exceptional_cases(self) -> tuple[GenericCase, ...]:
        return tuple(case for case in self.cases if case.exceptional)

    def to_sympy(self) -> sp.Expr:
        return self.generic_formula

    def sample_points(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return tuple(case.sample_point() for case in self.generic_cases)

    def as_cell_set(self) -> CellSet:
        if self.generic_split is not None and not self.parameters:
            return CellSet(
                self.variables,
                self.generic_split.generic_cells,
                self.cad_result.cad.cells_by_level,
                self.generic_formula,
            )
        cells: list[CADCell] = []
        for case in self.generic_cases:
            cells.extend(case.solution_cells)
        all_vars = (*self.parameters, *self.variables)
        return CellSet(
            all_vars, tuple(cells), self.cad_result.cad.cells_by_level, self.generic_formula
        )

    def as_function(self) -> GenericCADFunction:
        return GenericCADFunction(self)


@dataclass(frozen=True)
class GenericCADFunction:
    """Evaluable generic decomposition over parameter space."""

    result: GenericCADResult

    def __call__(
        self, point: Mapping[sp.Symbol, sp.Expr | int | float] | Sequence[sp.Expr | int | float]
    ) -> sp.Expr:
        case = self.case_for(point)
        if case is None or not case.has_solution:
            return sp.false
        return simplify_qe_formula(
            case.solution_formula.subs(self._subs(point)), implication_minimize=False
        )

    def case_for(
        self, point: Mapping[sp.Symbol, sp.Expr | int | float] | Sequence[sp.Expr | int | float]
    ) -> GenericCase | None:
        subs = self._subs(point)
        for case in self.result.cases:
            value = sp.simplify(case.param_condition.subs(subs))
            if value == sp.true or value is sp.true:
                return case
            try:
                if bool(value):
                    return case
            except TypeError:
                pass
        return None

    def exceptional(
        self, point: Mapping[sp.Symbol, sp.Expr | int | float] | Sequence[sp.Expr | int | float]
    ) -> bool:
        case = self.case_for(point)
        return bool(case and case.exceptional)

    def _subs(
        self, point: Mapping[sp.Symbol, sp.Expr | int | float] | Sequence[sp.Expr | int | float]
    ) -> dict[sp.Symbol, sp.Expr | int | float]:
        params = self.result.parameters
        if isinstance(point, Mapping):
            subs = {param: point[param] for param in params if param in point}
        else:
            if len(point) != len(params):
                raise ValueError(
                    f"expected {len(params)} parameter coordinate(s), got {len(point)}"
                )
            subs = dict(zip(params, point, strict=True))
        if len(subs) != len(params):
            missing = [sp.sstr(param) for param in params if param not in subs]
            raise ValueError(f"missing parameter coordinate(s): {missing}")
        return subs


def _sample_mapping(
    symbols: Sequence[sp.Symbol], cell: CADCell | None
) -> Mapping[sp.Symbol, sp.Expr]:
    if cell is None:
        return {}
    return {sym: sample_to_expr(sample) for sym, sample in zip(symbols, cell.sample, strict=True)}


def _cell_condition(
    cell: CADCell, variables: Sequence[sp.Symbol], cells_by_level: Mapping[int, Sequence[CADCell]]
) -> sp.Expr:
    return cells_to_formula((cell,), tuple(variables), cells_by_level)


def _is_full_dim(cell: CADCell, cells_by_level: Mapping[int, Sequence[CADCell]]) -> bool:
    return _cell_dim_with_levels(cell, cells_by_level) == cell.level


def _cells_above_param(
    cells: Sequence[CADCell], param_cell: CADCell | None, param_count: int
) -> tuple[CADCell, ...]:
    if param_count == 0:
        return tuple(cells)
    assert param_cell is not None
    prefix = param_cell.index
    return tuple(cell for cell in cells if cell.index[:param_count] == prefix)


def _case_solution_formula(
    cells: Sequence[CADCell],
    all_symbols: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
) -> sp.Expr:
    if not cells:
        return sp.false
    return simplify_qe_formula(
        cells_to_formula(cells, tuple(all_symbols), cells_by_level), implication_minimize=False
    )


def _make_param_cases(
    cad_result: CADResult,
    parameters: tuple[sp.Symbol, ...],
    variables: tuple[sp.Symbol, ...],
    base_formula: sp.Expr,
) -> tuple[GenericCase, ...]:
    cad_obj = cad_result.cad
    selected = tuple(cad_result.cells)
    param_count = len(parameters)
    if param_count == 0:
        solution = base_formula
        return (
            GenericCase(
                param_condition=sp.true,
                solution_formula=solution,
                param_cell=None,
                solution_cells=selected,
                param_sample={},
                dimension=0,
                exceptional=False,
            ),
        )
    cases: list[GenericCase] = []
    for param_cell in cad_obj.cells_by_level.get(param_count, tuple()):
        cells = _cells_above_param(selected, param_cell, param_count)
        param_condition = simplify_qe_formula(
            _cell_condition(param_cell, parameters, cad_obj.cells_by_level),
            implication_minimize=False,
        )
        cell_formula = _case_solution_formula(
            cells, (*parameters, *variables), cad_obj.cells_by_level
        )
        solution = (
            sp.false
            if not cells
            else simplify_qe_formula(
                sp.And(param_condition, base_formula, cell_formula, evaluate=False),
                implication_minimize=False,
            )
        )
        full_dim = _is_full_dim(param_cell, cad_obj.cells_by_level)
        cases.append(
            GenericCase(
                param_condition=param_condition,
                solution_formula=solution,
                param_cell=param_cell,
                solution_cells=cells,
                param_sample=_sample_mapping(parameters, param_cell),
                dimension=_cell_dim_with_levels(param_cell, cad_obj.cells_by_level),
                exceptional=not full_dim,
            )
        )
    return tuple(cases)


def _generic_formula(cases: Sequence[GenericCase]) -> sp.Expr:
    return _formula_from_conditions(
        case.param_condition
        for case in cases
        if not case.exceptional and case.solution_formula != sp.false
    )


def _exceptional_formula(cases: Sequence[GenericCase]) -> sp.Expr:
    return _formula_from_conditions(case.param_condition for case in cases if case.exceptional)


def generic_cad(
    formula: sp.Expr | Formula,
    variables: Sequence[sp.Symbol | str],
    *,
    parameters: Sequence[sp.Symbol | str] = (),
    output: GenericOutput = "formula",
    strategy: str = "auto",
    domain: str = "reals",
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    strict: bool = False,
    return_result: bool = True,
):
    """Compute a generic cylindrical decomposition over parameter space.

    Parameters are placed before ordinary variables in the underlying CAD. The
    generic formula is the union of full-dimensional parameter cells with
    nonempty selected fibers. Lower-dimensional parameter cells are reported as
    exceptional cases, including exceptional cells with empty fibers.
    """

    dom = normalize_domain(domain)
    if dom.value != "reals":
        if strict:
            raise NotImplementedError("generic CAD currently supports only the real domain")
        var_tuple = _normalize_variables(variables)
        param_tuple = _normalize_variables(parameters)
        result = GenericCADResult(
            generic_formula=sp.false,
            exceptional_formula=sp.false,
            cases=(),
            variables=var_tuple,
            parameters=param_tuple,
            cad_result=cad(sp.false, (*param_tuple, *var_tuple), return_result=True),
            output=output,
            status="unknown",
            diagnostics={"reason": f"unsupported generic CAD domain {dom.value}"},
            generic_split=None,
            exceptional_causes=(),
        )
        return result if return_result else result.generic_formula
    var_tuple = _normalize_variables(variables)
    param_tuple = _normalize_variables(parameters)
    all_symbols = (*param_tuple, *var_tuple)
    if isinstance(formula, Formula):
        base_raw = to_sympy(formula)
    else:
        base_raw = (
            formula
            if isinstance(formula, (sp.Basic, sp.logic.boolalg.Boolean))
            else sp.sympify(formula)
        )
    formula_with_assumptions = apply_assumptions(base_raw, assumptions)
    base_expr, _ = _normalize_formula(formula_with_assumptions)
    cad_result = cad(
        formula_with_assumptions,
        all_symbols,
        output="cells",
        strategy=strategy,
        domain=dom.value,
        assumptions=None,
    )
    cases = _make_param_cases(cad_result, param_tuple, var_tuple, base_expr)
    all_causes = (
        *input_boundary_causes(parse_formula(base_expr), all_symbols),
        *projection_causes(cad_result.cad, all_symbols),
    )
    generic_split = None
    generic_formula = simplify_qe_formula(_generic_formula(cases), implication_minimize=False)
    exceptional_formula = simplify_qe_formula(
        _exceptional_formula(cases), implication_minimize=False
    )
    exceptional_causes = tuple(cause for cause in all_causes if cause.source == "input_boundary")
    if not param_tuple:
        boundary_causes = tuple(cause for cause in all_causes if cause.source == "input_boundary")
        if (
            isinstance(
                base_expr, (sp.LessThan, sp.StrictLessThan, sp.GreaterThan, sp.StrictGreaterThan)
            )
            and boundary_causes
        ):
            delta = sp.expand(base_expr.lhs - base_expr.rhs)
            if isinstance(base_expr, (sp.LessThan, sp.StrictLessThan)):
                generic_formula = sp.StrictLessThan(delta, 0)
            else:
                generic_formula = sp.StrictGreaterThan(delta, 0)
            exceptional_formula = sp.Eq(delta, 0)
            generic_split = GenericSplit(
                generic_cells=tuple(cad_result.cells),
                exceptional_cells=tuple(cad_result.cells),
                exceptional_polys=tuple(cause.polynomial for cause in boundary_causes),
                exceptional_causes=boundary_causes,
                generic_formula=generic_formula,
                exceptional_formula=exceptional_formula,
            )
            exceptional_causes = boundary_causes
        else:
            generic_split = generic_split_from_cells(
                tuple(cad_result.cells), cad_result.cad, var_tuple, all_causes
            )
            generic_formula = generic_split.generic_formula
            exceptional_formula = generic_split.exceptional_formula
            exceptional_causes = generic_split.exceptional_causes
    result = GenericCADResult(
        generic_formula=simplify_qe_formula(generic_formula, implication_minimize=False),
        exceptional_formula=simplify_qe_formula(exceptional_formula, implication_minimize=False),
        cases=cases,
        variables=var_tuple,
        parameters=param_tuple,
        cad_result=cad_result,
        output=output,
        status="complete",
        diagnostics={
            "strategy": strategy,
            "domain": dom.value,
            "assumptions": tuple(map(sp.sstr, normalize_assumptions(assumptions))),
            "parameter_count": len(param_tuple),
            "variable_count": len(var_tuple),
            "generic_case_count": sum(not case.exceptional for case in cases),
            "generic_solution_case_count": sum(
                not case.exceptional and case.has_solution for case in cases
            ),
            "exceptional_case_count": sum(case.exceptional for case in cases),
            "exceptional_polynomial_count": len(exceptional_causes),
            "exceptional_causes": tuple(cause.source for cause in exceptional_causes),
            "cad_status": cad_result.status,
        },
        generic_split=generic_split,
        exceptional_causes=tuple(exceptional_causes),
    )
    if return_result:
        return result
    if output == "formula":
        return result.generic_formula
    if output == "cases":
        return result.cases
    if output == "cells":
        return result.as_cell_set()
    if output == "function":
        return result.as_function()
    raise ValueError(f"unsupported generic CAD output: {output!r}")


def generic_cad_text(
    text: str,
    *,
    variables: Sequence[sp.Symbol | str],
    parameters: Sequence[sp.Symbol | str] = (),
    symbols: Mapping[str, sp.Symbol] | None = None,
    output: GenericOutput = "formula",
    strategy: str = "auto",
    domain: str = "reals",
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    strict: bool = False,
    return_result: bool = True,
):
    local_symbols = dict(symbols or {})
    for sym_like in (*parameters, *variables):
        if isinstance(sym_like, str):
            local_symbols.setdefault(sym_like, sp.Symbol(sym_like, real=True))
        else:
            local_symbols.setdefault(sym_like.name, sym_like)
    expr, _ = parse_formula_text(text, symbols=local_symbols)
    return generic_cad(
        expr,
        variables,
        parameters=parameters,
        output=output,
        strategy=strategy,
        domain=domain,
        assumptions=assumptions,
        strict=strict,
        return_result=return_result,
    )


__all__ = [
    "GenericCADFunction",
    "GenericCADResult",
    "GenericCase",
    "GenericOutput",
    "generic_cad",
    "generic_cad_text",
]
