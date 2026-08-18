from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal

import sympy as sp

from ..algebraic.comparison import compare_samples
from ..algebraic.samples import sample_to_expr
from ..cad.decomposition import CompleteCAD, decomp_collins_complete
from ..cad.lifting.stack import CADCell
from ..domains import apply_assumptions, normalize_assumptions, normalize_domain
from ..formula import Formula, formula_polynomials, parse_formula, parse_formula_text, to_sympy
from ..preprocess.semialgebraicize import semialgebraicize
from ..qe.complete import CellUnion, cells_to_formula, evaluate_formula_on_cell, qe_by_complete_cad
from ..reconstruct.cylindrical import path_condition
from ..simplify.result import simplify_qe_formula
from ..topology.operations import apply_topological_operation

CADOutput = Literal["formula", "cells", "function", "tree"]
FormulaForm = Literal["nested", "dnf"]
TopoOp = Literal["closure", "interior", "boundary", "exterior", "components"]


@dataclass(frozen=True)
class CADOptions:
    """Options accepted by :func:`cad`."""

    output: CADOutput = "formula"
    operation: TopoOp | None = None
    strategy: str = "auto"
    domain: str = "reals"
    assumptions: tuple[sp.Expr, ...] = ()
    max_cells: int | None = None
    timeout: float | None = None
    diagnostics: bool = True
    strict: bool = False
    formula_form: FormulaForm = "nested"
    max_formula_terms: int = 512


@dataclass(frozen=True)
class CellSet:
    """A finite set of CAD cells selected by a semialgebraic condition."""

    variables: tuple[sp.Symbol, ...]
    cells: tuple[CADCell, ...]
    cells_by_level: Mapping[int, tuple[CADCell, ...]] = field(default_factory=dict)
    formula: sp.Expr = sp.false

    def __len__(self) -> int:
        return len(self.cells)

    def __iter__(self):
        return iter(self.cells)

    @property
    def is_empty(self) -> bool:
        return not self.cells

    @property
    def indices(self) -> tuple[tuple[int, ...], ...]:
        return tuple(cell.index for cell in self.cells)

    def sample_points(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return tuple(
            {
                var: sample_to_expr(sample)
                for var, sample in zip(self.variables, cell.sample, strict=True)
            }
            for cell in self.cells
        )

    def to_formula(self) -> sp.Expr:
        return self.formula


@dataclass(frozen=True)
class CADTreeNode:
    """Node in an evaluable cylindrical decomposition tree."""

    cell: CADCell
    truth: bool | None
    children: tuple[CADTreeNode, ...] = ()

    @property
    def level(self) -> int:
        return self.cell.level

    @property
    def index(self) -> tuple[int, ...]:
        return self.cell.index


@dataclass(frozen=True)
class CADFunction:
    """Compact evaluable representation of a cylindrical decomposition.

    The object stores the CAD tree and selected full-dimensional cells. Point
    membership first locates the CAD leaf containing the point and then returns
    the stored truth value for that leaf. This avoids repeatedly evaluating a
    reconstructed Boolean formula, while preserving the formula as a fallback
    display representation.
    """

    variables: tuple[sp.Symbol, ...]
    formula: sp.Expr
    cad: CompleteCAD
    cell_set: CellSet
    operation: str | None = None
    tree: tuple[CADTreeNode, ...] = ()

    def __post_init__(self) -> None:
        if not self.tree:
            object.__setattr__(self, "tree", build_cad_tree(self.cad, self.cell_set))

    def __call__(self, *values) -> bool:
        return self.contains(values)

    def contains(
        self, point: Sequence[sp.Expr | int | float] | Mapping[sp.Symbol, sp.Expr | int | float]
    ) -> bool:
        subs = self._point_subs(point)
        leaf = self.locate_cell(subs)
        if leaf is not None:
            selected = {cell.index for cell in self.cell_set.cells}
            if leaf.index in selected:
                return True
            node = self.tree_by_index().get(leaf.index)
            if node is not None and node.truth is not None:
                return node.truth
        value = self.formula.subs(subs)
        return bool(value)

    def locate_cell(
        self,
        point: Sequence[sp.Expr | int | float] | Mapping[sp.Symbol, sp.Expr | int | float],
    ) -> CADCell | None:
        subs = self._point_subs(point)
        parent_index: tuple[int, ...] | None = None
        found: CADCell | None = None
        for level, var in enumerate(self.variables, start=1):
            value = sp.sympify(subs[var])
            candidates = [
                cell
                for cell in self.cad.cells_by_level.get(level, ())
                if cell.parent_index == parent_index
            ]
            found = _cell_containing_value(candidates, value)
            if found is None:
                return None
            parent_index = found.index
        return found

    def sample_points(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return self.cell_set.sample_points()

    def to_formula(self) -> sp.Expr:
        return self.formula

    def restrict(self, assignments: Mapping[sp.Symbol, sp.Expr | int | float]) -> CADFunction:
        remaining = tuple(var for var in self.variables if var not in assignments)
        restricted = sp.simplify(self.formula.subs(assignments))
        if not remaining:
            restricted_cell_set = CellSet((), tuple(), {}, restricted)
            return CADFunction((), restricted, self.cad, restricted_cell_set, self.operation)
        result = cad(restricted, remaining, output="function")
        return result.as_function()

    def project_cell_set(self, variables: Sequence[sp.Symbol]) -> CellSet:
        level = len(tuple(variables))
        if tuple(variables) != self.variables[:level]:
            raise NotImplementedError(
                "projection currently requires a prefix of the decomposition variables"
            )
        projected: dict[tuple[int, ...], CADCell] = {}
        for cell in self.cell_set.cells:
            key = cell.index[:level]
            if key:
                match = next(c for c in self.cad.cells_by_level[level] if c.index == key)
                projected[key] = match
        cells = tuple(projected[key] for key in sorted(projected))
        formula = (
            cells_to_formula(cells, tuple(variables), self.cad.cells_by_level)
            if cells
            else sp.false
        )
        return CellSet(
            tuple(variables), cells, self.cad.cells_by_level, simplify_qe_formula(formula)
        )

    def project(self, variables: Sequence[sp.Symbol]) -> CADFunction:
        """Return a compact CAD function for a prefix projection.

        The current implementation projects only onto a prefix of the CAD
        variable order. The projected formula is rebuilt from selected prefix
        cells and then decomposed in the lower-dimensional space so subsequent
        membership queries use a compact function object rather than a raw cell
        set.
        """

        cell_set = self.project_cell_set(variables)
        return cad(cell_set.formula, tuple(variables), output="function").as_function()

    def tree_by_index(self) -> Mapping[tuple[int, ...], CADTreeNode]:
        nodes: dict[tuple[int, ...], CADTreeNode] = {}

        def visit(node: CADTreeNode) -> None:
            nodes[node.index] = node
            for child in node.children:
                visit(child)

        for root in self.tree:
            visit(root)
        return nodes

    def _point_subs(
        self, point: Sequence[sp.Expr | int | float] | Mapping[sp.Symbol, sp.Expr | int | float]
    ) -> dict[sp.Symbol, sp.Expr | int | float]:
        if isinstance(point, Mapping):
            subs = {var: point[var] for var in self.variables if var in point}
        else:
            if len(point) != len(self.variables):
                raise ValueError(f"expected {len(self.variables)} coordinate(s), got {len(point)}")
            subs = dict(zip(self.variables, point, strict=True))
        if len(subs) != len(self.variables):
            missing = [sp.sstr(var) for var in self.variables if var not in subs]
            raise ValueError(f"missing coordinate(s): {missing}")
        return subs


@dataclass(frozen=True)
class CADResult:
    """Public result for CAD requests."""

    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    cad: CompleteCAD
    cell_set: CellSet
    output: CADOutput
    operation: str | None
    status: str
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    function: CADFunction | None = None

    @property
    def cells(self) -> tuple[CADCell, ...]:
        return self.cell_set.cells

    @property
    def tree(self) -> tuple[CADTreeNode, ...]:
        return self.as_function().tree

    def as_cell_set(self) -> CellSet:
        return self.cell_set

    def cell_count_by_level(self) -> Mapping[int, int]:
        return self.cad.cell_count_by_level()

    def proj_poly_count_by_level(self) -> Mapping[int, int]:
        return self.cad.proj_poly_count_by_level()

    def sample_points(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return self.cell_set.sample_points()

    def to_sympy(self) -> sp.Expr:
        return self.formula

    def to_dnf(self) -> sp.Expr:
        return sp.simplify_logic(self.formula, form="dnf")

    def as_function(self) -> CADFunction:
        return self.function or CADFunction(
            self.variables, self.formula, self.cad, self.cell_set, self.operation
        )


def _normalize_variables(variables: Sequence[sp.Symbol | str]) -> tuple[sp.Symbol, ...]:
    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for var in variables:
        sym = sp.Symbol(var, real=True) if isinstance(var, str) else var
        if sym not in seen:
            out.append(sym)
            seen.add(sym)
    return tuple(out)


def _normalize_formula(expr_or_formula: sp.Expr | Formula) -> tuple[sp.Expr, Formula]:
    if isinstance(expr_or_formula, (sp.Basic, sp.logic.boolalg.Boolean)):
        expr = expr_or_formula
        return expr, parse_formula(expr)
    return to_sympy(expr_or_formula), expr_or_formula


def _full_cells(cad_obj: CompleteCAD, variables: Sequence[sp.Symbol]) -> tuple[CADCell, ...]:
    return cad_obj.cells_by_level.get(len(variables), tuple())


def _cell_dim_with_levels(cell: CADCell, cells_by_level: Mapping[int, Sequence[CADCell]]) -> int:
    dim = 0
    for level in range(1, cell.level + 1):
        prefix = cell.index[:level]
        ancestor = next(c for c in cells_by_level[level] if c.index == prefix)
        if ancestor.kind == "sector":
            dim += 1
    return dim


def _closed_cell_condition(
    cell: CADCell, variables: Sequence[sp.Symbol], cells_by_level: Mapping[int, Sequence[CADCell]]
) -> sp.Expr:
    return path_condition(cell, variables, cells_by_level, closed=True)


def _open_cell_condition(
    cell: CADCell, variables: Sequence[sp.Symbol], cells_by_level: Mapping[int, Sequence[CADCell]]
) -> sp.Expr:
    if any(
        next(c for c in cells_by_level[level] if c.index == cell.index[:level]).kind == "section"
        for level in range(1, cell.level + 1)
    ):
        return sp.false
    return path_condition(cell, variables, cells_by_level, closed=False)


def _formula_from_conditions(conditions: Iterable[sp.Expr]) -> sp.Expr:
    kept = [cond for cond in conditions if cond is not sp.false and cond != sp.false]
    if not kept:
        return sp.false
    return simplify_qe_formula(sp.simplify_logic(sp.Or(*kept), form="dnf"))


def _operation_formula(
    selected_cells: Sequence[CADCell],
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    operation: str | None,
) -> sp.Expr:
    if operation is None:
        return cells_to_formula(selected_cells, variables, cells_by_level)
    if operation == "closure":
        return _formula_from_conditions(
            _closed_cell_condition(cell, variables, cells_by_level) for cell in selected_cells
        )
    if operation == "interior":
        full_dim = len(variables)
        return _formula_from_conditions(
            _open_cell_condition(cell, variables, cells_by_level)
            for cell in selected_cells
            if _cell_dim_with_levels(cell, cells_by_level) == full_dim
        )
    if operation == "boundary":
        closure = _operation_formula(selected_cells, variables, cells_by_level, "closure")
        interior = _operation_formula(selected_cells, variables, cells_by_level, "interior")
        return simplify_qe_formula(sp.simplify_logic(sp.And(closure, sp.Not(interior)), form="dnf"))
    raise ValueError(f"unsupported topological operation: {operation!r}")


def _select_formula_cells(
    cad_obj: CompleteCAD, formula: Formula, variables: Sequence[sp.Symbol]
) -> tuple[CADCell, ...]:
    return tuple(
        cell
        for cell in _full_cells(cad_obj, variables)
        if evaluate_formula_on_cell(formula, cell, variables)
    )


def _build_cad_for_formula(formula: Formula, variables: Sequence[sp.Symbol]) -> CompleteCAD:
    polys = formula_polynomials(formula)
    if not polys:
        polys = [sp.Integer(1)]
    return decomp_collins_complete(polys, variables)


def _cell_containing_value(cells: Sequence[CADCell], value: sp.Expr) -> CADCell | None:
    for cell in cells:
        left, right = cell.interval or (None, None)
        if cell.kind == "section":
            if left is not None and sp.simplify(value - sample_to_expr(left)) == 0:
                return cell
            continue
        lower_ok = True if left is None else _value_gt_sample(value, left)
        upper_ok = True if right is None else _value_lt_sample(value, right)
        if lower_ok and upper_ok:
            return cell
    return None


def _value_gt_sample(value: sp.Expr, sample) -> bool:
    try:
        if hasattr(sample, "to_expr") or hasattr(sample, "value") or hasattr(sample, "expr"):
            return compare_samples(_expr_sample(value), sample) > 0
    except Exception:
        pass
    return bool(sp.N(value - sample_to_expr(sample)) > 0)


def _value_lt_sample(value: sp.Expr, sample) -> bool:
    try:
        if hasattr(sample, "to_expr") or hasattr(sample, "value") or hasattr(sample, "expr"):
            return compare_samples(_expr_sample(value), sample) < 0
    except Exception:
        pass
    return bool(sp.N(value - sample_to_expr(sample)) < 0)


def _expr_sample(value: sp.Expr):
    from ..algebraic.samples import RationalSample

    rational = sp.Rational(value) if value.is_Rational or value.is_Integer else None
    if rational is None:
        raise TypeError("only rational point-location values use exact sample comparison")
    return RationalSample(rational)


def build_cad_tree(cad_obj: CompleteCAD, cell_set: CellSet) -> tuple[CADTreeNode, ...]:
    selected = {cell.index for cell in cell_set.cells}
    by_parent: dict[tuple[int, ...] | None, list[CADCell]] = {}
    for level_cells in cad_obj.cells_by_level.values():
        for cell in level_cells:
            by_parent.setdefault(cell.parent_index, []).append(cell)

    def build(cell: CADCell) -> CADTreeNode:
        children = tuple(
            build(child)
            for child in sorted(by_parent.get(cell.index, ()), key=lambda c: c.stack_position)
        )
        truth = cell.index in selected if not children else None
        return CADTreeNode(cell=cell, truth=truth, children=children)

    roots = sorted(by_parent.get(None, ()), key=lambda c: c.stack_position)
    return tuple(build(root) for root in roots)


def _make_result(
    expr: sp.Expr,
    formula: Formula,
    variables: tuple[sp.Symbol, ...],
    options: CADOptions,
) -> CADResult:
    cad_obj = _build_cad_for_formula(formula, variables)
    base_cells = _select_formula_cells(cad_obj, formula, variables)
    if options.operation == "components":
        topo_cells = tuple(sorted(base_cells, key=lambda cell: cell.index))
        raw_formula = cells_to_formula(
            topo_cells,
            variables,
            cad_obj.cells_by_level,
            form=options.formula_form,
            max_terms=options.max_formula_terms,
        )
    elif options.operation is None:
        topo_cells = tuple(sorted(base_cells, key=lambda cell: cell.index))
        raw_formula = cells_to_formula(
            topo_cells,
            variables,
            cad_obj.cells_by_level,
            form=options.formula_form,
            max_terms=options.max_formula_terms,
        )
    else:
        topo = apply_topological_operation(base_cells, cad_obj, variables, options.operation)
        topo_cells = topo.cells
        raw_formula = topo.formula
    public_formula = simplify_qe_formula(raw_formula, implication_minimize=False)
    output_formula = public_formula
    if options.operation is None:
        cell_union = CellUnion(
            variables=variables,
            cells=base_cells,
            formula=raw_formula,
            cells_by_level=cad_obj.cells_by_level,
        )
        output_formula = simplify_qe_formula(
            raw_formula, cell_union=cell_union, implication_minimize=False
        )
    selected_cells = topo_cells
    cell_set = CellSet(variables, selected_cells, cad_obj.cells_by_level, output_formula)
    function = CADFunction(variables, output_formula, cad_obj, cell_set, options.operation)
    diagnostics = {
        "strategy": options.strategy,
        "domain": options.domain,
        "cell_count_by_level": cad_obj.cell_count_by_level(),
        "projection_polynomial_count_by_level": cad_obj.proj_poly_count_by_level(),
        "operation": options.operation,
        "assumptions": tuple(map(sp.sstr, options.assumptions)),
        "input_formula": expr,
        "formula_form": options.formula_form,
        "max_formula_terms": options.max_formula_terms,
    }
    return CADResult(
        formula=output_formula,
        variables=variables,
        cad=cad_obj,
        cell_set=cell_set,
        output=options.output,
        operation=options.operation,
        status="complete",
        diagnostics=diagnostics,
        function=function,
    )


def cad(
    formula: sp.Expr | Formula,
    variables: Sequence[sp.Symbol | str],
    *,
    output: CADOutput = "formula",
    operation: TopoOp | None = None,
    strategy: str = "auto",
    domain: str = "reals",
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    max_cells: int | None = None,
    timeout: float | None = None,
    diagnostics: bool = True,
    strict: bool = False,
    return_result: bool = True,
    formula_form: FormulaForm = "nested",
    max_formula_terms: int = 512,
):
    """Compute a cylindrical algebraic decomposition for a real formula."""

    dom = normalize_domain(domain)
    variables = _normalize_variables(variables)
    if dom.value != "reals":
        if strict:
            raise NotImplementedError("CAD currently supports only the real domain")
        empty_cad = decomp_collins_complete([sp.Integer(1)], variables)
        empty_set = CellSet(variables, tuple(), empty_cad.cells_by_level, sp.false)
        result = CADResult(
            sp.false,
            variables,
            empty_cad,
            empty_set,
            output,
            operation,
            "unknown",
            {"reason": f"unsupported CAD domain {dom.value}"},
        )
        return result if return_result else result.formula
    base_expr = (
        to_sympy(formula)
        if not isinstance(formula, (sp.Basic, sp.logic.boolalg.Boolean))
        else formula
    )
    expr, normalized = _normalize_formula(apply_assumptions(base_expr, assumptions))
    options = CADOptions(
        output=output,
        operation=operation,
        strategy=strategy,
        domain=dom.value,
        assumptions=normalize_assumptions(assumptions),
        max_cells=max_cells,
        timeout=timeout,
        diagnostics=diagnostics,
        strict=strict,
        formula_form=formula_form,
        max_formula_terms=max_formula_terms,
    )
    prep = semialgebraicize(expr, variables=variables)
    if prep.changed and prep.aux_vars:
        if len(prep.aux_vars) > 3:
            result = _make_result(sp.true, parse_formula(sp.true), variables, options)
            qe_formula = sp.true
        else:
            internal_vars = tuple(variables) + tuple(prep.aux_vars)
            qe_result = qe_by_complete_cad(
                internal_vars,
                tuple(("exists", aux) for aux in prep.aux_vars),
                prep.formula,
                free_variables=variables,
            )
            qe_formula = qe_result.formula
            result = _make_result(qe_formula, parse_formula(qe_formula), variables, options)
        diag = dict(result.diagnostics)
        diag.update(
            {
                "preprocessed": True,
                "preprocess_aux_vars": tuple(map(sp.sstr, prep.aux_vars)),
                "preprocess_notes": prep.notes,
                "preprocessed_formula": prep.sympy_expr,
                "qe_formula": qe_formula,
                "preprocess_elimination_limited": len(prep.aux_vars) > 3,
            }
        )
        result = CADResult(
            result.formula,
            result.variables,
            result.cad,
            result.cell_set,
            result.output,
            result.operation,
            result.status,
            diag,
            result.function,
        )
    elif prep.changed:
        result = _make_result(prep.sympy_expr, prep.formula, variables, options)
    else:
        result = _make_result(expr, normalized, variables, options)
    if return_result:
        return result
    if output == "formula":
        return result.formula
    if output == "cells":
        return result.cell_set
    if output == "function":
        return result.as_function()
    if output == "tree":
        return result.tree
    raise ValueError(f"unsupported CAD output: {output!r}")


def cad_text(
    text: str,
    *,
    variables: Sequence[sp.Symbol | str] | None = None,
    symbols: Mapping[str, sp.Symbol] | None = None,
    output: CADOutput = "formula",
    operation: TopoOp | None = None,
    strategy: str = "auto",
    domain: str = "reals",
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    max_cells: int | None = None,
    timeout: float | None = None,
    diagnostics: bool = True,
    strict: bool = False,
    return_result: bool = True,
    formula_form: FormulaForm = "nested",
    max_formula_terms: int = 512,
):
    local_symbols = dict(symbols or {})
    if variables is not None:
        for var in variables:
            if isinstance(var, str):
                local_symbols.setdefault(var, sp.Symbol(var, real=True))
            else:
                local_symbols.setdefault(var.name, var)
    expr, formula = parse_formula_text(text, symbols=local_symbols)
    if variables is None:
        variables = tuple(sorted(expr.free_symbols, key=lambda s: s.name))
    return cad(
        formula,
        variables,
        output=output,
        operation=operation,
        strategy=strategy,
        domain=domain,
        assumptions=assumptions,
        max_cells=max_cells,
        timeout=timeout,
        diagnostics=diagnostics,
        strict=strict,
        return_result=return_result,
        formula_form=formula_form,
        max_formula_terms=max_formula_terms,
    )


__all__ = [
    "CADFunction",
    "CADOptions",
    "CADOutput",
    "FormulaForm",
    "CADResult",
    "CADTreeNode",
    "CellSet",
    "TopoOp",
    "build_cad_tree",
    "cad",
    "cad_text",
]
