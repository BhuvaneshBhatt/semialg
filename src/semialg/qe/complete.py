from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..cad.decomposition import CompleteCAD, decomp_collins_complete
from ..cad.lifting.stack import CADCell
from ..formula import And, Atom, BoolConst, Formula, Not, Or, formula_polynomials, to_sympy
from ..reconstruct.merge import compressed_formula_from_cells, dnf_formula_from_cells
from ..simplify.result import simplify_qe_formula
from .blocks import QuantifierBlock, quantifiers_to_blocks


@dataclass(frozen=True)
class CellUnion:
    """A finite union of cells in a CAD projection space.

    The formula field is a conservative reconstruction of the same union. The
    cell list remains the semantic object; callers can simplify or pretty-print
    the formula later without losing the audited CAD provenance.
    """

    variables: tuple[sp.Symbol, ...]
    cells: tuple[CADCell, ...]
    formula: sp.Expr
    cells_by_level: Mapping[int, tuple[CADCell, ...]] = field(default_factory=dict)

    @property
    def is_empty(self) -> bool:
        return not self.cells

    @property
    def cell_indices(self) -> tuple[tuple[int, ...], ...]:
        return tuple(cell.index for cell in self.cells)


@dataclass(frozen=True)
class QEDiagnostics:
    """Small audit record for the complete-CAD QE driver."""

    requested_variables: tuple[sp.Symbol, ...]
    internal_variables: tuple[sp.Symbol, ...]
    free_variables: tuple[sp.Symbol, ...]
    quantified_variables: tuple[sp.Symbol, ...]
    quantifier_blocks: tuple[QuantifierBlock, ...]
    variable_reordered: bool
    projection_level: int
    full_cell_count: int
    projected_cell_count: int
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class CompleteQEResult:
    """Result returned by the conservative Collins-based QE path.

    The result intentionally keeps both a reconstructed formula and the CAD
    objects used to derive it. ``formula`` is the public quantifier-free result;
    ``cell_union`` is the semantic projection object for non-sentence results;
    ``cell_truth`` records the final truth value on each free-variable cell.
    """

    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    free_variables: tuple[sp.Symbol, ...]
    quantified_variables: tuple[sp.Symbol, ...]
    quantifiers: tuple[tuple[str, sp.Symbol], ...]
    cad: CompleteCAD
    quantifier_blocks: tuple[QuantifierBlock, ...] = ()
    cell_union: CellUnion | None = None
    backend: str = "collins-complete-qe"
    status: str = "complete"
    is_sentence: bool = False
    truth_value: bool | None = None
    satisfying_cell_indices: tuple[tuple[int, ...], ...] = ()
    cell_truth: Mapping[tuple[int, ...], bool] = field(default_factory=dict)
    witness_samples: tuple[Mapping[str, sp.Expr], ...] = ()
    diagnostics: QEDiagnostics | None = None


def _sign_key(expr: sp.Expr) -> str:
    return sp.sstr(sp.expand(expr))


def norm_quant_map(quantifiers: Sequence[tuple[str, sp.Symbol]]) -> dict[sp.Symbol, str]:
    out: dict[sp.Symbol, str] = {}
    for qname, sym in quantifiers:
        q = qname.lower()
        if q not in {"exists", "forall"}:
            raise ValueError(f"unsupported quantifier: {qname!r}")
        if sym in out:
            raise ValueError(f"variable {sym} is quantified more than once")
        out[sym] = q
    return out


def _ordered_unique(symbols: Sequence[sp.Symbol]) -> tuple[sp.Symbol, ...]:
    seen: set[sp.Symbol] = set()
    out: list[sp.Symbol] = []
    for sym in symbols:
        if sym not in seen:
            out.append(sym)
            seen.add(sym)
    return tuple(out)


def _matrix_symbols(matrix: Formula) -> tuple[sp.Symbol, ...]:
    return tuple(sorted(to_sympy(matrix).free_symbols, key=lambda s: s.name))


def norm_internal_order(
    vars_: Sequence[sp.Symbol],
    quantifiers: Sequence[tuple[str, sp.Symbol]],
    matrix: Formula,
    free_variables: Sequence[sp.Symbol] | None = None,
) -> tuple[tuple[sp.Symbol, ...], tuple[sp.Symbol, ...], tuple[sp.Symbol, ...], tuple[str, ...]]:
    """Return ``(internal, free, quantified, notes)`` for complete-CAD QE.

    Collins projection eliminates suffix variables by propagating truth values
    upward through the CAD tree. Public callers may provide variables in any
    order, so the complete driver normalizes to:

        free variables, then quantified variables in prefix order.

    This preserves quantifier semantics while removing the old user-facing
    order restriction.
    """

    qmap = norm_quant_map(quantifiers)
    quantified = tuple(sym for _, sym in quantifiers)
    matrix_syms = _matrix_symbols(matrix)
    requested = _ordered_unique(tuple(vars_) + matrix_syms + quantified)
    qset = set(quantified)

    if free_variables is None:
        free = tuple(sym for sym in requested if sym not in qset)
    else:
        free = _ordered_unique(tuple(free_variables))
        unknown = [sym for sym in free if sym not in requested]
        if unknown:
            raise ValueError(f"free variable(s) not present in formula/variable list: {unknown!r}")
        if any(sym in qset for sym in free):
            raise ValueError("a variable cannot be both free and quantified")

    missing_matrix = [sym for sym in matrix_syms if sym not in set(free) | qset]
    if missing_matrix:
        free = free + tuple(sym for sym in missing_matrix if sym not in free)

    internal = free + quantified
    notes: list[str] = []
    if tuple(vars_) != internal:
        notes.append(
            "internal variable order normalized to free variables followed by quantified prefix order"
        )
    # Touch qmap so duplicate/invalid quantifier checks are guaranteed even when
    # quantified is empty; the returned ordering itself uses the prefix sequence.
    _ = qmap
    return internal, free, quantified, tuple(notes)


def _atom_truth(
    atom: Atom, signs: Mapping[str, int], sample: Sequence[object], variables: Sequence[sp.Symbol]
) -> bool:
    key = _sign_key(atom.expr)
    sign = signs.get(key)
    if sign is None:
        substitutions = {
            sym: sample_to_expr(val) for sym, val in zip(variables, sample, strict=True)
        }
        value = sp.expand(atom.expr).subs(substitutions)
        sign = int(sp.sign(value))
    if atom.op == "=":
        return sign == 0
    if atom.op == "!=":
        return sign != 0
    if atom.op == "<":
        return sign < 0
    if atom.op == "<=":
        return sign <= 0
    if atom.op == ">":
        return sign > 0
    if atom.op == ">=":
        return sign >= 0
    raise ValueError(f"unknown atomic operator: {atom.op!r}")


def evaluate_formula_on_cell(
    formula: Formula, cell: CADCell, variables: Sequence[sp.Symbol] | None = None
) -> bool:
    """Evaluate a quantifier-free formula using recorded signs on ``cell``."""

    if variables is None:
        variables = tuple(sorted(to_sympy(formula).free_symbols, key=lambda s: s.name))
    if isinstance(formula, BoolConst):
        return formula.value
    if isinstance(formula, Atom):
        return _atom_truth(formula, cell.signs, cell.sample, variables)
    if isinstance(formula, And):
        return all(evaluate_formula_on_cell(arg, cell, variables) for arg in formula.args)
    if isinstance(formula, Or):
        return any(evaluate_formula_on_cell(arg, cell, variables) for arg in formula.args)
    if isinstance(formula, Not):
        return not evaluate_formula_on_cell(formula.arg, cell, variables)
    raise TypeError(f"unsupported formula node: {type(formula)!r}")


def _combine(values: Sequence[bool], quantifier: str) -> bool:
    if quantifier == "exists":
        return any(values)
    if quantifier == "forall":
        return all(values)
    raise ValueError(f"unsupported quantifier: {quantifier!r}")


def _cell_condition(
    cell: CADCell, variables: Sequence[sp.Symbol], cells_by_level: Mapping[int, Sequence[CADCell]]
) -> sp.Expr:
    pieces: list[sp.Expr] = []
    for level in range(1, cell.level + 1):
        prefix = cell.index[:level]
        level_cell = next(c for c in cells_by_level[level] if c.index == prefix)
        var = variables[level - 1]
        left, right = level_cell.interval or (None, None)
        if left is not None and right is not None and left == right:
            pieces.append(sp.Eq(var, sample_to_expr(left)))
        else:
            if left is not None:
                pieces.append(var > sample_to_expr(left))
            if right is not None:
                pieces.append(var < sample_to_expr(right))
    if not pieces:
        return sp.true
    return sp.And(*pieces)


def cells_to_formula(
    cells: Sequence[CADCell],
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    form: str = "nested",
    max_terms: int = 512,
) -> sp.Expr:
    """Convert a union of projection cells into a conservative SymPy formula.

    ``form="nested"`` shares common CAD prefixes and merges contiguous stack
    blocks. ``form="dnf"`` preserves the older one-disjunct-per-cell output.
    """

    if not cells:
        return sp.false
    if form == "dnf":
        return sp.simplify_logic(
            dnf_formula_from_cells(cells, variables, cells_by_level), form="dnf"
        )
    result = compressed_formula_from_cells(cells, variables, cells_by_level, max_terms=max_terms)
    if result.stats.fallback_used:
        return sp.simplify_logic(result.formula, form="dnf")
    return result.formula


def _dummy_cad() -> CompleteCAD:
    return decomp_collins_complete([sp.Integer(1)], [sp.Symbol("_dummy", real=True)])


def _leaf_witnesses(
    cad: CompleteCAD,
    variables: Sequence[sp.Symbol],
    matrix: Formula,
    *,
    limit: int = 8,
) -> tuple[Mapping[str, sp.Expr], ...]:
    witnesses: list[Mapping[str, sp.Expr]] = []
    full_level = len(variables)
    if full_level == 0:
        return tuple()
    for cell in cad.cells_by_level.get(full_level, tuple()):
        if evaluate_formula_on_cell(matrix, cell, variables):
            witnesses.append(
                {
                    sp.sstr(sym): sample_to_expr(sample)
                    for sym, sample in zip(variables, cell.sample, strict=True)
                }
            )
            if len(witnesses) >= limit:
                break
    return tuple(witnesses)


def qe_by_complete_cad(
    vars_: Sequence[sp.Symbol],
    quantifiers: Sequence[tuple[str, sp.Symbol]],
    matrix: Formula,
    *,
    free_variables: Sequence[sp.Symbol] | None = None,
) -> CompleteQEResult:
    """Eliminate quantified real variables using the conservative CAD path.

    Public callers may provide variables in any order. The implementation uses
    an internal order with free variables first and quantified variables in their
    prefix order, which is the order needed for suffix elimination by CAD truth
    propagation. Arbitrary alternating quantifier prefixes are supported.
    """

    requested_variables = tuple(vars_)
    quantifiers = tuple((q.lower(), sym) for q, sym in quantifiers)
    qmap = norm_quant_map(quantifiers)
    variables, free, quantified, notes = norm_internal_order(
        requested_variables, quantifiers, matrix, free_variables
    )
    blocks = quantifiers_to_blocks(quantifiers)

    if not variables:
        truth = bool(to_sympy(matrix))
        diag = QEDiagnostics(
            requested_variables=requested_variables,
            internal_variables=variables,
            free_variables=free,
            quantified_variables=quantified,
            quantifier_blocks=blocks,
            variable_reordered=tuple(requested_variables) != variables,
            projection_level=0,
            full_cell_count=0,
            projected_cell_count=0,
            notes=notes,
        )
        return CompleteQEResult(
            formula=sp.true if truth else sp.false,
            variables=variables,
            free_variables=free,
            quantified_variables=quantified,
            quantifiers=quantifiers,
            quantifier_blocks=blocks,
            cad=_dummy_cad(),
            is_sentence=True,
            truth_value=truth,
            diagnostics=diag,
        )

    polys = formula_polynomials(matrix)
    if not polys:
        formula = to_sympy(matrix)
        truth = bool(formula) if quantified or not free else None
        cad = decomp_collins_complete([sp.Integer(1)], variables)
        cell_union = None
        if free:
            cells = cad.cells_by_level[len(free)]
            raw_cell_formula = formula if bool(formula) else sp.false
            cell_union = CellUnion(
                variables=free,
                cells=cells if bool(formula) else tuple(),
                formula=raw_cell_formula,
                cells_by_level={len(free): cells},
            )
            raw_cell_formula = simplify_qe_formula(raw_cell_formula, cell_union=cell_union)
            cell_union = CellUnion(
                variables=free,
                cells=cell_union.cells,
                formula=raw_cell_formula,
                cells_by_level=cell_union.cells_by_level,
            )
        diag = QEDiagnostics(
            requested_variables=requested_variables,
            internal_variables=variables,
            free_variables=free,
            quantified_variables=quantified,
            quantifier_blocks=blocks,
            variable_reordered=tuple(requested_variables) != variables,
            projection_level=len(free),
            full_cell_count=len(cad.cells),
            projected_cell_count=0 if cell_union is None else len(cell_union.cells),
            notes=notes,
        )
        return CompleteQEResult(
            formula=sp.true if truth else sp.false if truth is False else formula,
            variables=variables,
            free_variables=free,
            quantified_variables=quantified,
            quantifiers=quantifiers,
            quantifier_blocks=blocks,
            cad=cad,
            cell_union=cell_union,
            is_sentence=not free,
            truth_value=truth,
            diagnostics=diag,
        )

    cad = decomp_collins_complete(polys, variables)
    n = len(variables)
    free_level = len(free)
    current: dict[tuple[int, ...], bool] = {
        cell.index: evaluate_formula_on_cell(matrix, cell, variables)
        for cell in cad.cells_by_level[n]
    }

    for level in range(n, free_level, -1):
        var = variables[level - 1]
        quantifier = qmap[var]
        grouped: dict[tuple[int, ...], list[bool]] = {}
        for idx, value in current.items():
            grouped.setdefault(idx[:-1], []).append(value)
        current = {parent: _combine(values, quantifier) for parent, values in grouped.items()}

    witnesses = _leaf_witnesses(cad, variables, matrix)

    if free_level == 0:
        truth = current.get(tuple(), False)
        diag = QEDiagnostics(
            requested_variables=requested_variables,
            internal_variables=variables,
            free_variables=free,
            quantified_variables=quantified,
            quantifier_blocks=blocks,
            variable_reordered=tuple(requested_variables) != variables,
            projection_level=0,
            full_cell_count=len(cad.cells),
            projected_cell_count=1 if truth else 0,
            notes=notes,
        )
        return CompleteQEResult(
            formula=sp.true if truth else sp.false,
            variables=variables,
            free_variables=free,
            quantified_variables=quantified,
            quantifiers=quantifiers,
            quantifier_blocks=blocks,
            cad=cad,
            is_sentence=True,
            truth_value=truth,
            satisfying_cell_indices=(tuple(),) if truth else (),
            cell_truth=current,
            witness_samples=witnesses,
            diagnostics=diag,
        )

    free_cells = tuple(
        cell for cell in cad.cells_by_level[free_level] if current.get(cell.index, False)
    )
    raw_formula = cells_to_formula(free_cells, free, cad.cells_by_level)
    cell_union = CellUnion(
        variables=free,
        cells=free_cells,
        formula=raw_formula,
        cells_by_level={
            level: tuple(cells)
            for level, cells in cad.cells_by_level.items()
            if level <= free_level
        },
    )
    formula = simplify_qe_formula(raw_formula, cell_union=cell_union)
    cell_union = CellUnion(
        variables=free, cells=free_cells, formula=formula, cells_by_level=cell_union.cells_by_level
    )
    diag = QEDiagnostics(
        requested_variables=requested_variables,
        internal_variables=variables,
        free_variables=free,
        quantified_variables=quantified,
        quantifier_blocks=blocks,
        variable_reordered=tuple(requested_variables) != variables,
        projection_level=free_level,
        full_cell_count=len(cad.cells),
        projected_cell_count=len(free_cells),
        notes=notes,
    )
    return CompleteQEResult(
        formula=formula,
        variables=variables,
        free_variables=free,
        quantified_variables=quantified,
        quantifiers=quantifiers,
        quantifier_blocks=blocks,
        cad=cad,
        cell_union=cell_union,
        is_sentence=False,
        truth_value=None,
        satisfying_cell_indices=tuple(cell.index for cell in free_cells),
        cell_truth=current,
        witness_samples=witnesses,
        diagnostics=diag,
    )


__all__ = [
    "CellUnion",
    "CompleteQEResult",
    "QEDiagnostics",
    "cells_to_formula",
    "evaluate_formula_on_cell",
    "qe_by_complete_cad",
    "qe_from_cad",
]


def qe_from_cad(
    cad: CompleteCAD,
    vars_: Sequence[sp.Symbol],
    quantifiers: Sequence[tuple[str, sp.Symbol]],
    matrix: Formula,
    *,
    free_variables: Sequence[sp.Symbol] | None = None,
    backend: str | None = None,
) -> CompleteQEResult:
    """Run the QE truth-propagation layer over an already lifted CAD.

    This is used by the strategy planner when it selects an active reduced CAD
    backend. The caller is responsible for ensuring that ``cad`` is complete
    for the requested invariant; safe reduced backends do so by certification
    and fallback.
    """

    requested_variables = tuple(vars_)
    quantifiers = tuple((q.lower(), sym) for q, sym in quantifiers)
    qmap = norm_quant_map(quantifiers)
    variables, free, quantified, notes = norm_internal_order(
        requested_variables, quantifiers, matrix, free_variables
    )
    if tuple(cad.tower.variables) != tuple(variables):
        raise ValueError("CAD variable order does not match normalized QE order")
    blocks = quantifiers_to_blocks(quantifiers)
    n = len(variables)
    free_level = len(free)
    if n == 0:
        truth = bool(to_sympy(matrix))
        return CompleteQEResult(
            formula=sp.true if truth else sp.false,
            variables=variables,
            free_variables=free,
            quantified_variables=quantified,
            quantifiers=quantifiers,
            quantifier_blocks=blocks,
            cad=cad,
            backend=backend or cad.backend,
            is_sentence=True,
            truth_value=truth,
        )
    current: dict[tuple[int, ...], bool] = {
        cell.index: evaluate_formula_on_cell(matrix, cell, variables)
        for cell in cad.cells_by_level[n]
    }
    for level in range(n, free_level, -1):
        var = variables[level - 1]
        quantifier = qmap[var]
        grouped: dict[tuple[int, ...], list[bool]] = {}
        for idx, value in current.items():
            grouped.setdefault(idx[:-1], []).append(value)
        current = {parent: _combine(values, quantifier) for parent, values in grouped.items()}
    witnesses = _leaf_witnesses(cad, variables, matrix)
    if free_level == 0:
        truth = current.get(tuple(), False)
        diag = QEDiagnostics(
            requested_variables=requested_variables,
            internal_variables=variables,
            free_variables=free,
            quantified_variables=quantified,
            quantifier_blocks=blocks,
            variable_reordered=tuple(requested_variables) != variables,
            projection_level=0,
            full_cell_count=len(cad.cells),
            projected_cell_count=1 if truth else 0,
            notes=notes,
        )
        return CompleteQEResult(
            formula=sp.true if truth else sp.false,
            variables=variables,
            free_variables=free,
            quantified_variables=quantified,
            quantifiers=quantifiers,
            quantifier_blocks=blocks,
            cad=cad,
            backend=backend or cad.backend,
            is_sentence=True,
            truth_value=truth,
            satisfying_cell_indices=(tuple(),) if truth else (),
            cell_truth=current,
            witness_samples=witnesses,
            diagnostics=diag,
        )
    free_cells = tuple(
        cell for cell in cad.cells_by_level[free_level] if current.get(cell.index, False)
    )
    raw_formula = cells_to_formula(free_cells, free, cad.cells_by_level)
    cell_union = CellUnion(
        variables=free,
        cells=free_cells,
        formula=raw_formula,
        cells_by_level={
            level: tuple(cells)
            for level, cells in cad.cells_by_level.items()
            if level <= free_level
        },
    )
    formula = simplify_qe_formula(raw_formula, cell_union=cell_union)
    cell_union = CellUnion(
        variables=free, cells=free_cells, formula=formula, cells_by_level=cell_union.cells_by_level
    )
    diag = QEDiagnostics(
        requested_variables=requested_variables,
        internal_variables=variables,
        free_variables=free,
        quantified_variables=quantified,
        quantifier_blocks=blocks,
        variable_reordered=tuple(requested_variables) != variables,
        projection_level=free_level,
        full_cell_count=len(cad.cells),
        projected_cell_count=len(free_cells),
        notes=notes,
    )
    return CompleteQEResult(
        formula=formula,
        variables=variables,
        free_variables=free,
        quantified_variables=quantified,
        quantifiers=quantifiers,
        quantifier_blocks=blocks,
        cad=cad,
        backend=backend or cad.backend,
        cell_union=cell_union,
        is_sentence=False,
        truth_value=None,
        satisfying_cell_indices=tuple(cell.index for cell in free_cells),
        cell_truth=current,
        witness_samples=witnesses,
        diagnostics=diag,
    )
