from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..algebraic.signs import sign_at_sample
from ..cad_algorithms.decomposition import (
    CompleteCAD,
    decomp_collins_complete,
    try_decomp_groebner_variety,
)
from ..cad_algorithms.lifting.stack import CADCell
from ..cad_algorithms.projection.groebner import build_groebner_variety_projection
from ..context import with_computation_context
from ..exact_arithmetic import exact_truth
from ..formula import (
    And,
    Atom,
    BoolConst,
    Formula,
    Not,
    Or,
    formula_polynomials,
    parse_formula,
    to_sympy,
)
from ..heuristics import suggest_cad_variable_order
from ..presolve import fourier_motzkin_eliminate, presolve_semialgebraic
from ..simplify.formula import simplify_qe_formula
from ..structural_keys import symbol_identity_key
from .blocks import QuantifierBlock, quantifiers_to_blocks
from .cell_reconstruction import cells_to_formula, finite_variety_formula


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
    projection_polynomials: Mapping[int, tuple[sp.Expr, ...]] = field(default_factory=dict)

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
    witness_samples: tuple[Mapping[sp.Symbol, sp.Expr], ...] = ()
    diagnostics: QEDiagnostics | None = None


def _sign_key(expr: sp.Expr) -> str:
    return sp.srepr(sp.expand(expr))


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


def _canonicalize_same_kind_quantifier_blocks(
    quantifiers: Sequence[tuple[str, sp.Symbol]],
) -> tuple[tuple[str, sp.Symbol], ...]:
    """Canonicalize variable order inside contiguous equal-kind blocks.

    Variables within one existential or universal block may be permuted without
    changing first-order semantics.  A deterministic order lets projection/CAD
    caches be reused across callers that present the same block in different
    permutations, while never moving a variable across a quantifier boundary.
    """

    out: list[tuple[str, sp.Symbol]] = []
    start = 0
    quantifiers = tuple(quantifiers)
    while start < len(quantifiers):
        qname = quantifiers[start][0]
        end = start + 1
        while end < len(quantifiers) and quantifiers[end][0] == qname:
            end += 1
        block = sorted(quantifiers[start:end], key=lambda item: symbol_identity_key(item[1]))
        out.extend(block)
        start = end
    return tuple(out)


def _ordered_unique(symbols: Sequence[sp.Symbol]) -> tuple[sp.Symbol, ...]:
    seen: set[sp.Symbol] = set()
    out: list[sp.Symbol] = []
    for sym in symbols:
        if sym not in seen:
            out.append(sym)
            seen.add(sym)
    return tuple(out)


def _matrix_symbols(matrix: Formula) -> tuple[sp.Symbol, ...]:
    return tuple(sorted(to_sympy(matrix).free_symbols, key=symbol_identity_key))


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

    This preserves quantifier semantics while accepting caller variable orders
    independently of the internal projection order.
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
        try:
            poly = (
                sp.Poly(sp.expand(atom.expr), *variables, domain="EX")
                if variables
                else sp.Poly(sp.expand(atom.expr))
            )
            sign = sign_at_sample(poly, sample)
        except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
            substitutions = {
                sym: sample_to_expr(val) for sym, val in zip(variables, sample, strict=True)
            }
            value = sp.expand(atom.expr).subs(substitutions)
            exact = sp.sign(value)
            if exact not in (-1, 0, 1):
                raise ValueError(
                    f"could not determine exact atom sign for {sp.sstr(atom.expr)}"
                ) from None
            sign = int(exact)
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
        variables = tuple(sorted(to_sympy(formula).free_symbols, key=symbol_identity_key))
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


def _dummy_cad() -> CompleteCAD:
    return decomp_collins_complete([sp.Integer(1)], [sp.Symbol("_dummy", real=True)])


def _leaf_witnesses(
    cad: CompleteCAD,
    variables: Sequence[sp.Symbol],
    matrix: Formula,
    *,
    limit: int = 8,
) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
    witnesses: list[Mapping[sp.Symbol, sp.Expr]] = []
    full_level = len(variables)
    if full_level == 0:
        return tuple()
    for cell in cad.cells_by_level.get(full_level, tuple()):
        if evaluate_formula_on_cell(matrix, cell, variables):
            witnesses.append(
                {
                    sym: sample_to_expr(sample)
                    for sym, sample in zip(variables, cell.sample, strict=True)
                }
            )
            if len(witnesses) >= limit:
                break
    return tuple(witnesses)


def _projection_polynomials_for_free_levels(
    cad: CompleteCAD, free_level: int
) -> dict[int, tuple[sp.Expr, ...]]:
    return {
        level: tuple(poly.as_expr() for poly in cad.tower.level(level).polynomials)
        for level in range(1, free_level + 1)
    }


def _preferred_free_atoms(matrix: Formula, free: tuple[sp.Symbol, ...]) -> tuple[sp.Expr, ...]:
    expr = to_sympy(matrix)
    if not isinstance(expr, sp.And):
        return ()
    free_set = set(free)
    return tuple(
        atom
        for atom in expr.args
        if getattr(atom, "is_Relational", False) and atom.free_symbols <= free_set
    )


def _auto_qe_order(
    polys: Sequence[sp.Expr],
    free: tuple[sp.Symbol, ...],
    quantifiers: tuple[tuple[str, sp.Symbol], ...],
    *,
    strategy: str,
) -> tuple[tuple[sp.Symbol, ...], tuple[tuple[str, sp.Symbol], ...]]:
    """Reorder only semantically interchangeable QE variable blocks."""

    if strategy in {"preserve", "none"} or not polys:
        return free, quantifiers
    free_order = (
        suggest_cad_variable_order(polys, free, strategy=strategy).order if len(free) > 1 else free
    )
    reordered: list[tuple[str, sp.Symbol]] = []
    start = 0
    while start < len(quantifiers):
        qname = quantifiers[start][0]
        end = start + 1
        while end < len(quantifiers) and quantifiers[end][0] == qname:
            end += 1
        block = tuple(sym for _, sym in quantifiers[start:end])
        if len(block) > 1:
            block = suggest_cad_variable_order(polys, block, strategy=strategy).order
        reordered.extend((qname, sym) for sym in block)
        start = end
    return tuple(free_order), tuple(reordered)


@with_computation_context
def _qe_by_complete_cad_result(
    vars_: Sequence[sp.Symbol],
    quantifiers: Sequence[tuple[str, sp.Symbol]],
    matrix: Formula,
    *,
    free_variables: Sequence[sp.Symbol] | None = None,
    variable_order_strategy: str = "auto",
    use_presolve: bool = True,
    allow_variety_cad: bool = True,
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
    note_list = list(notes)

    if variable_order_strategy in {"preserve", "none"} and quantifiers:
        canonical_quantifiers = _canonicalize_same_kind_quantifier_blocks(quantifiers)
        if canonical_quantifiers != quantifiers:
            quantifiers = canonical_quantifiers
            quantified = tuple(sym for _, sym in quantifiers)
            variables = (*free, *quantified)
            qmap = norm_quant_map(quantifiers)
            note_list.append(
                "canonicalized variable order within equal-kind quantifier blocks for CAD reuse"
            )

    if use_presolve and variables:
        # Only the innermost existential block may be substituted directly in
        # the quantifier-free matrix.  Eliminating an outer existential through
        # an equality that depends on a later universal variable would change
        # ``exists x. forall y`` into the invalid ``forall y. exists x`` semantics.
        existential: tuple[sp.Symbol, ...] = ()
        if quantifiers and quantifiers[-1][0] == "exists":
            start = len(quantifiers) - 1
            while start > 0 and quantifiers[start - 1][0] == "exists":
                start -= 1
            existential = tuple(sym for _, sym in quantifiers[start:])
        presolved = presolve_semialgebraic(to_sympy(matrix), variables, eliminate=existential)
        if presolved.changed:
            removed = {var for var, _ in presolved.substitutions}
            matrix = parse_formula(presolved.formula)
            quantifiers = tuple((qname, sym) for qname, sym in quantifiers if sym not in removed)
            free = tuple(sym for sym in free if sym not in removed)
            quantified = tuple(sym for _, sym in quantifiers)
            variables = tuple(sym for sym in variables if sym not in removed)
            qmap = norm_quant_map(quantifiers)
            if removed:
                note_list.append(
                    "presolve eliminated affine existential variables: "
                    + ", ".join(sym.name for sym in removed)
                )

        # After affine substitution, a remaining innermost existential block
        # may be a purely linear conjunction.  Fourier-Motzkin removes it
        # exactly without constructing projection polynomials in those vars.
        trailing_exists: tuple[sp.Symbol, ...] = ()
        if quantifiers and quantifiers[-1][0] == "exists":
            start = len(quantifiers) - 1
            while start > 0 and quantifiers[start - 1][0] == "exists":
                start -= 1
            trailing_exists = tuple(sym for _, sym in quantifiers[start:])
        if trailing_exists:
            fm = fourier_motzkin_eliminate(to_sympy(matrix), trailing_exists)
            if fm is not None:
                fm_formula, fm_removed = fm
                matrix = parse_formula(fm_formula)
                removed_set = set(fm_removed)
                quantifiers = tuple(
                    (qname, sym) for qname, sym in quantifiers if sym not in removed_set
                )
                quantified = tuple(sym for _, sym in quantifiers)
                variables = tuple(sym for sym in variables if sym not in removed_set)
                qmap = norm_quant_map(quantifiers)
                note_list.append(
                    "Fourier-Motzkin eliminated linear existential variables: "
                    + ", ".join(sym.name for sym in fm_removed)
                )

    polys_for_order = formula_polynomials(matrix)
    if variable_order_strategy not in {"preserve", "none"} and variables:
        free, quantifiers = _auto_qe_order(
            polys_for_order, free, quantifiers, strategy=variable_order_strategy
        )
        quantified = tuple(sym for _, sym in quantifiers)
        variables = (*free, *quantified)
        qmap = norm_quant_map(quantifiers)
        if tuple(variables) != tuple(sym for sym in requested_variables if sym in set(variables)):
            note_list.append(
                f"CAD variable ordering selected by {variable_order_strategy!r} heuristic"
            )
    notes = tuple(note_list)
    blocks = quantifiers_to_blocks(quantifiers)
    preferred_atoms = _preferred_free_atoms(matrix, free)

    if not variables:
        truth = exact_truth(to_sympy(matrix))
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
                projection_polynomials=_projection_polynomials_for_free_levels(cad, len(free)),
            )
            raw_cell_formula = simplify_qe_formula(
                raw_cell_formula, cell_union=cell_union, preferred_atoms=preferred_atoms
            )
            cell_union = CellUnion(
                variables=free,
                cells=cell_union.cells,
                formula=raw_cell_formula,
                cells_by_level=cell_union.cells_by_level,
                projection_polynomials=cell_union.projection_polynomials,
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

    groebner_projection = None
    if allow_variety_cad and all(qname == "exists" for qname, _ in quantifiers):
        groebner_projection = build_groebner_variety_projection(matrix, variables)
    if groebner_projection is not None:
        variety_cad = try_decomp_groebner_variety(groebner_projection)
        if variety_cad is None:
            cad = decomp_collins_complete(polys, variables)
            note_list.append(
                "Groebner variety lifting could not certify a fiber; used complete Collins CAD"
            )
        else:
            cad = variety_cad
            note_list.append(
                "used zero-dimensional Groebner variety projection/lifting for existential CAD"
            )
        notes = tuple(note_list)
    else:
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
            backend=(
                "groebner-variety-qe"
                if cad.backend == "groebner-variety"
                else "collins-complete-qe"
            ),
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
    if cad.tower.metadata.get("variety_only", False):
        raw_formula = finite_variety_formula(free_cells, free)
    else:
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
        projection_polynomials=_projection_polynomials_for_free_levels(cad, free_level),
    )
    if cad.tower.metadata.get("variety_only", False):
        formula = simplify_qe_formula(raw_formula, preferred_atoms=preferred_atoms)
    else:
        formula = simplify_qe_formula(
            raw_formula, cell_union=cell_union, preferred_atoms=preferred_atoms
        )
    cell_union = CellUnion(
        variables=free,
        cells=free_cells,
        formula=formula,
        cells_by_level=cell_union.cells_by_level,
        projection_polynomials=cell_union.projection_polynomials,
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
        backend=(
            "groebner-variety-qe" if cad.backend == "groebner-variety" else "collins-complete-qe"
        ),
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


@with_computation_context
def qe_by_complete_cad(
    vars_: Sequence[sp.Symbol],
    quantifiers: Sequence[tuple[str, sp.Symbol]],
    matrix: Formula,
    *,
    free_variables: Sequence[sp.Symbol] | None = None,
    variable_order_strategy: str = "auto",
    use_presolve: bool = True,
    allow_variety_cad: bool = True,
    return_result: bool = False,
) -> sp.Expr | CompleteQEResult:
    """Eliminate quantified real variables and return the resulting formula.

    Set ``return_result=True`` to retain the complete CAD decomposition,
    truth-propagation metadata, witnesses, and diagnostics.
    """

    result = _qe_by_complete_cad_result(
        vars_,
        quantifiers,
        matrix,
        free_variables=free_variables,
        variable_order_strategy=variable_order_strategy,
        use_presolve=use_presolve,
        allow_variety_cad=allow_variety_cad,
    )
    return result if return_result else result.formula


def _qe_from_cad_result(
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
    preferred_atoms = _preferred_free_atoms(matrix, free)
    n = len(variables)
    free_level = len(free)
    if n == 0:
        truth = exact_truth(to_sympy(matrix))
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
        projection_polynomials=_projection_polynomials_for_free_levels(cad, free_level),
    )
    formula = simplify_qe_formula(
        raw_formula, cell_union=cell_union, preferred_atoms=preferred_atoms
    )
    cell_union = CellUnion(
        variables=free,
        cells=free_cells,
        formula=formula,
        cells_by_level=cell_union.cells_by_level,
        projection_polynomials=cell_union.projection_polynomials,
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


def qe_from_cad(
    cad: CompleteCAD,
    vars_: Sequence[sp.Symbol],
    quantifiers: Sequence[tuple[str, sp.Symbol]],
    matrix: Formula,
    *,
    free_variables: Sequence[sp.Symbol] | None = None,
    backend: str | None = None,
    return_result: bool = False,
) -> sp.Expr | CompleteQEResult:
    """Run QE truth propagation over an existing CAD and return the formula.

    Set ``return_result=True`` for CAD cells, witnesses, truth metadata, and
    diagnostics in a :class:`CompleteQEResult`.
    """

    result = _qe_from_cad_result(
        cad,
        vars_,
        quantifiers,
        matrix,
        free_variables=free_variables,
        backend=backend,
    )
    return result if return_result else result.formula
