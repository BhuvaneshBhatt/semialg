from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..cad_algorithms.lifting.stack import CADCell
from ..context import with_computation_context
from ..domains import apply_assumptions, normalize_assumptions, normalize_domain
from ..formula import Formula, parse_formula, parse_formula_text, to_sympy
from ..qe.complete import cells_to_formula
from ..simplify.formula import simplify_qe_formula
from ..symbol_resolution import build_symbol_table
from ._parametric_support import (
    ParametricCADBoundaryAnalysis,
    ParametricCADBoundaryCause,
    ParametricCADSplit,
    analyze_parametric_boundaries,
    input_boundary_causes,
    parametric_split_from_cells,
    projection_causes,
)
from .cylindrical import (
    CADResult,
    CellSet,
    _cell_dim_with_levels,
    _formula_from_conditions,
    _normalize_formula,
    _normalize_variables,
    cad,
)

ParametricCADOutput = Literal["result", "formula", "cases", "cells", "function"]


@dataclass(frozen=True)
class ParametricCADCase:
    """One parameter cell in a parametric CAD analysis."""

    param_condition: sp.Expr
    index: int
    solution_formula: sp.Expr
    param_cell: CADCell | None
    solution_cells: tuple[CADCell, ...]
    param_sample: Mapping[sp.Symbol, sp.Expr]
    dimension: int
    exceptional: bool = False
    variables: tuple[sp.Symbol, ...] = ()
    specialize_fibers: bool = True

    @property
    def has_solution(self) -> bool:
        return self.solution_cells != () and self.solution_formula != sp.false

    @property
    def condition(self) -> sp.Expr:
        return self.param_condition

    @property
    def sample(self) -> Mapping[sp.Symbol, sp.Expr]:
        return self.param_sample

    @property
    def parameter_cell(self) -> CADCell | None:
        return self.param_cell

    @property
    def specialized_formula(self) -> sp.Expr:
        if not self.has_solution:
            return sp.false
        return simplify_qe_formula(
            self.solution_formula.subs(dict(self.param_sample)), implication_minimize=False
        )

    @property
    def solution(self):
        """Representative cylindrical fiber, when exact extraction succeeds."""
        if not self.has_solution or not self.specialize_fibers:
            return None
        try:
            from ..cad_algorithms.cells import extract_cylindrical_solution

            return extract_cylindrical_solution(
                self.specialized_formula, self.variables, selected_only=True
            )
        except (NotImplementedError, ValueError, TypeError):
            return None

    def sample_point(self) -> Mapping[sp.Symbol, sp.Expr]:
        return self.param_sample

    def fiber_formula(self) -> sp.Expr:
        if not self.has_solution:
            return sp.false
        return simplify_qe_formula(
            self.solution_formula.subs(dict(self.param_sample)), implication_minimize=False
        )


@dataclass(frozen=True)
class ParametricCADResult:
    """Result returned by :func:`parametric_cad`."""

    generic_formula: sp.Expr
    exceptional_formula: sp.Expr
    cases: tuple[ParametricCADCase, ...]
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    cad_result: CADResult
    output: ParametricCADOutput
    status: str
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    generic_split: ParametricCADSplit | None = None
    exceptional_causes: tuple[ParametricCADBoundaryCause, ...] = ()
    exceptional_analysis: ParametricCADBoundaryAnalysis | None = None

    @property
    def strata(self) -> tuple[ParametricCADCase, ...]:
        return self.cases

    @property
    def stratum_count(self) -> int:
        return len(self.cases)

    @property
    def parameter_domain(self) -> sp.Expr:
        return _formula_from_conditions(case.param_condition for case in self.cases)

    @property
    def parameter_condition(self) -> sp.Expr:
        """Parameter values for which the original system has a nonempty fiber."""
        return _formula_from_conditions(
            case.param_condition for case in self.cases if case.has_solution
        )

    @property
    def generic_condition(self) -> sp.Expr:
        return self.generic_formula

    @property
    def exceptional_condition(self) -> sp.Expr:
        return self.exceptional_formula

    @property
    def generic_cases(self) -> tuple[ParametricCADCase, ...]:
        return tuple(case for case in self.cases if not case.exceptional and case.has_solution)

    @property
    def all_generic_cases(self) -> tuple[ParametricCADCase, ...]:
        return tuple(case for case in self.cases if not case.exceptional)

    @property
    def exceptional_cases(self) -> tuple[ParametricCADCase, ...]:
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

    def as_function(self) -> ParametricCADFunction:
        return ParametricCADFunction(self)

    def as_stratified_result(self):
        """Return parameter strata as guarded certified values."""
        from ..conditional import ConditionalBranch, conditional_result

        branches = [ConditionalBranch(case.param_condition, case) for case in self.cases]
        return conditional_result(
            self.parameters,
            branches,
            coverage_condition=self.parameter_domain,
            complete=self.status == "complete",
            disjoint=True,
            certified=self.status == "complete",
            method="parametric_cad",
            diagnostics={"stratum_count": len(self.cases)},
            normalize=False,
        )


@dataclass(frozen=True)
class ParametricCADFunction:
    """Evaluable generic decomposition over parameter space."""

    result: ParametricCADResult

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
    ) -> ParametricCADCase | None:
        subs = self._subs(point)
        for case in self.result.cases:
            value = sp.simplify(case.param_condition.subs(subs))
            if value == sp.true or value is sp.true:
                return case
            if value == sp.false or value is sp.false:
                continue
            if value.free_symbols:
                raise ValueError(
                    "parameter specialization did not produce an exact point condition"
                )
            # Do not use Python truth coercion for exact algebraic relations.
            # Reuse semialg's certified real decision layer for the closed
            # condition when SymPy has not simplified it to True/False.
            from ..decision import is_satisfiable

            if is_satisfiable(value, (), strategy="cad"):
                return case
        return None

    def specialize(
        self, point: Mapping[sp.Symbol, sp.Expr | int | float] | Sequence[sp.Expr | int | float]
    ) -> sp.Expr:
        """Return the exact fiber formula at one parameter point.

        The stored CAD is reused; specialization never rebuilds a decomposition.
        """
        return self(point)

    def normal_form(self) -> sp.Expr:
        """Return the exact full formula represented by the reusable CAD cases."""
        pieces = tuple(case.solution_formula for case in self.result.cases if case.has_solution)
        if not pieces:
            return sp.false
        return simplify_qe_formula(sp.Or(*pieces), implication_minimize=False)

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
    if param_cell is None:
        raise ValueError("parameter cell is required when param_count is nonzero")
    prefix = param_cell.index
    return tuple(cell for cell in cells if cell.index[:param_count] == prefix)


def _case_solution_formula(
    cells: Sequence[CADCell],
    all_symbols: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
) -> sp.Expr:
    """Reconstruct a parameterized fiber formula without rewriting root binders.

    Parametric CAD sections may contain ``root_of(p(a, x), x, k)`` expressions.
    Algebraic simplifiers that solve another conjunct for ``a`` and substitute
    it into the root-defining polynomial can collapse ``p`` to zero and destroy
    the ordered-root identity.  The CAD reconstruction is already an exact
    formula, so keep it structurally intact here; callers may simplify scalar
    conjuncts after parameter specialization.
    """
    if not cells:
        return sp.false
    return cells_to_formula(cells, tuple(all_symbols), cells_by_level)


def _make_param_cases(
    cad_result: CADResult,
    parameters: tuple[sp.Symbol, ...],
    variables: tuple[sp.Symbol, ...],
    base_formula: sp.Expr,
    *,
    specialize_fibers: bool = True,
) -> tuple[ParametricCADCase, ...]:
    """Construct generic and exceptional parameter cases from projection conditions."""
    cad_obj = cad_result.cad
    selected = tuple(cad_result.cells)
    param_count = len(parameters)
    if param_count == 0:
        solution = base_formula
        return (
            ParametricCADCase(
                param_condition=sp.true,
                index=0,
                solution_formula=solution,
                param_cell=None,
                solution_cells=selected,
                param_sample={},
                dimension=0,
                exceptional=False,
                variables=variables,
                specialize_fibers=specialize_fibers,
            ),
        )
    cases: list[ParametricCADCase] = []
    for case_index, param_cell in enumerate(cad_obj.cells_by_level.get(param_count, tuple())):
        cells = _cells_above_param(selected, param_cell, param_count)
        param_condition = simplify_qe_formula(
            _cell_condition(param_cell, parameters, cad_obj.cells_by_level),
            implication_minimize=False,
        )
        cell_formula = _case_solution_formula(
            cells, (*parameters, *variables), cad_obj.cells_by_level
        )
        # ``cell_formula`` already reconstructs exactly the selected CAD fiber
        # over this parameter cell, including the parameter-cell condition.
        # Re-simplifying it together with the original equation can substitute
        # solved parameter expressions into a bound ``root_of`` selector and
        # collapse its defining polynomial to zero.
        solution = sp.false if not cells else cell_formula
        full_dim = _is_full_dim(param_cell, cad_obj.cells_by_level)
        cases.append(
            ParametricCADCase(
                param_condition=param_condition,
                index=case_index,
                solution_formula=solution,
                param_cell=param_cell,
                solution_cells=cells,
                param_sample=_sample_mapping(parameters, param_cell),
                dimension=_cell_dim_with_levels(param_cell, cad_obj.cells_by_level),
                exceptional=not full_dim,
                variables=variables,
                specialize_fibers=specialize_fibers,
            )
        )
    return tuple(cases)


def _generic_formula(cases: Sequence[ParametricCADCase]) -> sp.Expr:
    return _formula_from_conditions(
        case.param_condition
        for case in cases
        if not case.exceptional and case.solution_formula != sp.false
    )


def _exceptional_formula(cases: Sequence[ParametricCADCase]) -> sp.Expr:
    return _formula_from_conditions(case.param_condition for case in cases if case.exceptional)


@with_computation_context
def parametric_cad(
    formula: sp.Expr | Formula,
    variables: Sequence[sp.Symbol | str],
    *,
    parameters: Sequence[sp.Symbol | str] = (),
    output: ParametricCADOutput = "result",
    strategy: str = "auto",
    domain: str = "reals",
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    strict: bool = False,
    specialize_fibers: bool = True,
):
    """Compute a parametric cylindrical decomposition over parameter space.

    The default return is the complete :class:`ParametricCADResult`. Use an
    explicit ``output`` selector only for a projected convenience view.

    Parameters are placed before ordinary variables in the underlying CAD. The
    generic formula is the union of full-dimensional parameter cells with
    nonempty selected fibers. Lower-dimensional parameter cells are reported as
    exceptional cases, including exceptional cells with empty fibers.
    """

    dom = normalize_domain(domain)
    if dom.value != "reals":
        if strict:
            raise NotImplementedError("parametric CAD supports only the real domain")
        var_tuple = _normalize_variables(variables)
        param_tuple = _normalize_variables(parameters)
        result = ParametricCADResult(
            generic_formula=sp.false,
            exceptional_formula=sp.false,
            cases=(),
            variables=var_tuple,
            parameters=param_tuple,
            cad_result=cad(sp.false, (*param_tuple, *var_tuple), return_result=True),
            output=output,
            status="unknown",
            diagnostics={"reason": f"unsupported parametric CAD domain {dom.value}"},
            generic_split=None,
            exceptional_causes=(),
        )
        if output == "result":
            return result
        raise NotImplementedError(
            f"parametric CAD supports only the real domain, not {dom.value!r}"
        )
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
        return_result=True,
    )
    cases = _make_param_cases(
        cad_result, param_tuple, var_tuple, base_expr, specialize_fibers=specialize_fibers
    )
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
            generic_split = ParametricCADSplit(
                generic_cells=tuple(cad_result.cells),
                exceptional_cells=tuple(cad_result.cells),
                exceptional_polys=tuple(cause.polynomial for cause in boundary_causes),
                exceptional_causes=boundary_causes,
                generic_formula=generic_formula,
                exceptional_formula=exceptional_formula,
            )
            exceptional_causes = boundary_causes
        else:
            generic_split = parametric_split_from_cells(
                tuple(cad_result.cells), cad_result.cad, var_tuple, all_causes
            )
            generic_formula = generic_split.generic_formula
            exceptional_formula = generic_split.exceptional_formula
            exceptional_causes = generic_split.exceptional_causes
    exceptional_analysis = None
    if param_tuple:
        exceptional_analysis = analyze_parametric_boundaries(base_expr, var_tuple, param_tuple)
    result = ParametricCADResult(
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
            "resource_limits": cad_result.diagnostics.get("resource_limits", {}),
            "planner": cad_result.diagnostics.get("planner", {}),
        },
        generic_split=generic_split,
        exceptional_causes=tuple(exceptional_causes),
        exceptional_analysis=exceptional_analysis,
    )
    if output == "result":
        return result
    if output == "formula":
        return result.generic_formula
    if output == "cases":
        return result.cases
    if output == "cells":
        return result.as_cell_set()
    if output == "function":
        return result.as_function()
    raise ValueError(f"unsupported parametric CAD output: {output!r}")


@with_computation_context
def parametric_cad_text(
    text: str,
    *,
    variables: Sequence[sp.Symbol | str],
    parameters: Sequence[sp.Symbol | str] = (),
    symbols: Mapping[str, sp.Symbol] | None = None,
    output: ParametricCADOutput = "result",
    strategy: str = "auto",
    domain: str = "reals",
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    strict: bool = False,
    specialize_fibers: bool = True,
):
    """Build a parametric cylindrical decomposition from a textual formula."""
    local_symbols = build_symbol_table(symbols, (*parameters, *variables))
    expr, _ = parse_formula_text(text, symbols=local_symbols)
    return parametric_cad(
        expr,
        variables,
        parameters=parameters,
        output=output,
        strategy=strategy,
        domain=domain,
        assumptions=assumptions,
        strict=strict,
        specialize_fibers=specialize_fibers,
    )


__all__ = [
    "ParametricCADFunction",
    "ParametricCADResult",
    "ParametricCADCase",
    "ParametricCADOutput",
    "parametric_cad",
    "parametric_cad_text",
]
