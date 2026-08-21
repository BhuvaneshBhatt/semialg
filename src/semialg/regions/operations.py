from __future__ import annotations

from collections.abc import Iterable, Sequence
from functools import cmp_to_key

import sympy as sp
from sympy.core.sympify import SympifyError
from sympy.logic.boolalg import Boolean
from sympy.polys.polyerrors import PolynomialError

from ..exact_arithmetic import compare_exact_reals
from ..implicit_geometry import decompose_cylindrical_formula_to_vertical_bounds_2d
from ..normalization import normalize_formula as _normalize_formula
from ..normalization import normalize_variables as _normalize_variables

FormulaLike = sp.Expr | Boolean | bool
_EXPECTED_ERRORS = (
    TypeError,
    ValueError,
    ArithmeticError,
    NotImplementedError,
    SympifyError,
    PolynomialError,
)


def _simplify_region(expr: sp.Expr) -> sp.Expr:
    if isinstance(expr, sp.core.relational.Relational):
        return expr
    try:
        if isinstance(expr, (sp.And, sp.Or, sp.Not)):
            return sp.simplify_logic(expr, form="dnf")
        return sp.simplify(expr)
    except _EXPECTED_ERRORS:
        return expr


def region_union(*regions: FormulaLike | Iterable[FormulaLike]) -> sp.Expr:
    """Return the union of implicit semialgebraic regions."""
    pieces = [_normalize_formula(region) for region in regions]
    return _simplify_region(sp.Or(*pieces) if pieces else sp.false)


def region_intersection(*regions: FormulaLike | Iterable[FormulaLike]) -> sp.Expr:
    """Return the intersection of implicit semialgebraic regions."""
    pieces = [_normalize_formula(region) for region in regions]
    return _simplify_region(sp.And(*pieces) if pieces else sp.true)


def region_complement(region: FormulaLike | Iterable[FormulaLike]) -> sp.Expr:
    """Return the complement of an implicit semialgebraic region."""
    return _simplify_region(sp.Not(_normalize_formula(region)))


def region_difference(
    lhs: FormulaLike | Iterable[FormulaLike],
    rhs: FormulaLike | Iterable[FormulaLike],
) -> sp.Expr:
    """Return ``lhs`` minus ``rhs`` for implicit semialgebraic regions."""
    return _simplify_region(sp.And(_normalize_formula(lhs), sp.Not(_normalize_formula(rhs))))


def _is_relational(expr: sp.Expr) -> bool:
    return isinstance(expr, sp.core.relational.Relational)


def _closure_atom(atom: sp.Expr) -> sp.Expr:
    if isinstance(atom, sp.StrictGreaterThan):
        return atom.lhs >= atom.rhs
    if isinstance(atom, sp.StrictLessThan):
        return atom.lhs <= atom.rhs
    return atom


def _interior_atom(atom: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr:
    if isinstance(atom, sp.GreaterThan):
        return atom.lhs > atom.rhs
    if isinstance(atom, sp.LessThan):
        return atom.lhs < atom.rhs
    if isinstance(atom, sp.Equality):
        diff = sp.simplify(atom.lhs - atom.rhs)
        return sp.true if diff == 0 else sp.false
    if isinstance(atom, sp.Unequality):
        diff = sp.simplify(atom.lhs - atom.rhs)
        return sp.false if diff == 0 else sp.true
    return atom


def _map_boolean(expr: sp.Expr, atom_fn) -> sp.Expr:
    if expr is sp.true or expr == sp.true:
        return sp.true
    if expr is sp.false or expr == sp.false:
        return sp.false
    if isinstance(expr, sp.And):
        return _simplify_region(sp.And(*[_map_boolean(arg, atom_fn) for arg in expr.args]))
    if isinstance(expr, sp.Or):
        return _simplify_region(sp.Or(*[_map_boolean(arg, atom_fn) for arg in expr.args]))
    if isinstance(expr, sp.Not):
        return _simplify_region(sp.Not(_map_boolean(expr.args[0], atom_fn)))
    if _is_relational(expr):
        return atom_fn(expr)
    return expr


def _syntactically_empty(expr: sp.Expr) -> bool:
    simplified = _simplify_region(expr)
    return simplified is sp.false or simplified == sp.false


def region_closure(
    region: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> sp.Expr:
    """Return a formula for the closure in common implicit cases.

    This first public wrapper relaxes strict inequalities to weak inequalities.
    That is exact for ordinary conjunctions/disjunctions of polynomial
    inequalities whose boundary is obtained by replacing strict comparisons by
    non-strict comparisons. Future improvements can replace this with CAD-cell
    adjacency closure for arbitrary Boolean formulas.
    """
    expr = _normalize_formula(region)
    if _syntactically_empty(expr):
        return sp.false
    return _simplify_region(_map_boolean(expr, _closure_atom))


def region_interior(
    region: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> sp.Expr:
    """Return a formula for the ambient interior in common implicit cases."""
    expr = _normalize_formula(region)
    vars_ = _normalize_variables(variables, expr)
    if _syntactically_empty(expr):
        return sp.false
    return _simplify_region(_map_boolean(expr, lambda atom: _interior_atom(atom, vars_)))


def _boundary_from_conjunction(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    atoms = list(expr.args) if isinstance(expr, sp.And) else [expr]
    if not atoms or not all(_is_relational(atom) for atom in atoms):
        return None
    boundaries: list[sp.Expr] = []
    for atom in atoms:
        if isinstance(atom, (sp.StrictLessThan, sp.StrictGreaterThan, sp.LessThan, sp.GreaterThan)):
            boundary_atom = sp.Eq(sp.simplify(atom.lhs - atom.rhs), 0)
            others = [region_closure(other, variables) for other in atoms if other is not atom]
            boundaries.append(sp.And(boundary_atom, *others))
    if not boundaries:
        return None
    return _simplify_region(sp.Or(*boundaries))


def _boundary_from_vertical_bounds_2d(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> sp.Expr | None:
    """Boundary formula from supported cylindrical 2D vertical-bound cells."""

    if len(variables) != 2:
        return None
    x, y = variables
    try:
        cells = decompose_cylindrical_formula_to_vertical_bounds_2d(expr, (x, y))
    except _EXPECTED_ERRORS:
        return None
    if not cells:
        return sp.false
    pieces: list[sp.Expr] = []
    for cell in cells:
        lo, hi = cell.x_interval
        if lo == -sp.oo and hi == sp.oo:
            x_open = sp.true
        elif lo == -sp.oo:
            x_open = x < hi
        elif hi == sp.oo:
            x_open = x > lo
        elif lo == hi:
            x_open = sp.Eq(x, lo)
        else:
            x_open = sp.And(x > lo, x < hi)
        for lower, upper in cell.y_bounds:
            if lower == upper:
                pieces.append(sp.And(sp.Eq(y, lower), x_open))
            else:
                if lower != -sp.oo:
                    pieces.append(sp.And(sp.Eq(y, lower), x_open))
                if upper != sp.oo:
                    pieces.append(sp.And(sp.Eq(y, upper), x_open))
                if lo != -sp.oo:
                    pieces.append(
                        sp.And(sp.Eq(x, lo), y >= lower.subs(x, lo), y <= upper.subs(x, lo))
                    )
                if hi != sp.oo and hi != lo:
                    pieces.append(
                        sp.And(sp.Eq(x, hi), y >= lower.subs(x, hi), y <= upper.subs(x, hi))
                    )
    return _simplify_region(sp.Or(*pieces) if pieces else sp.false)


def _boundary_from_structured_cad_cells_2d(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> sp.Expr | None:
    """Boundary formula from the complete CAD cell structure for 2D regions."""

    if len(variables) != 2:
        return None
    x, y = variables
    try:
        from ..cad.cells import extract_vertical_bounds_from_cad_2d

        cells = extract_vertical_bounds_from_cad_2d(expr, (x, y), full_dimensional_only=True)
    except _EXPECTED_ERRORS:
        return None
    if not cells:
        return None
    pieces: list[sp.Expr] = []
    for cell in cells:
        lo, hi = cell.x_interval
        if lo == -sp.oo and hi == sp.oo:
            x_open = sp.true
        elif lo == -sp.oo:
            x_open = x < hi
        elif hi == sp.oo:
            x_open = x > lo
        elif lo == hi:
            x_open = sp.Eq(x, lo)
        else:
            x_open = sp.And(x > lo, x < hi)
        for lower, upper in cell.y_bounds:
            if lower == upper:
                pieces.append(sp.And(sp.Eq(y, lower), x_open))
                continue
            if lower != -sp.oo:
                pieces.append(sp.And(sp.Eq(y, lower), x_open))
            if upper != sp.oo:
                pieces.append(sp.And(sp.Eq(y, upper), x_open))
            if lo != -sp.oo and lo != hi:
                pieces.append(sp.And(sp.Eq(x, lo), y >= lower.subs(x, lo), y <= upper.subs(x, lo)))
            if hi != sp.oo and hi != lo:
                pieces.append(sp.And(sp.Eq(x, hi), y >= lower.subs(x, hi), y <= upper.subs(x, hi)))
    return _simplify_region(sp.Or(*pieces) if pieces else sp.false)


def _boundary_from_cylindrical_solution_nd(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> sp.Expr | None:
    """Boundary candidate from arbitrary-dimensional cylindrical solution cells.

    For each full-dimensional cell, each finite coordinate bound contributes a
    face. Bounds at level k may depend on earlier variables, so faces are
    assembled from closed constraints at all levels with the selected coordinate
    replaced by an equality to the chosen bound. This is a conservative exact
    formula for the cell-boundary union when the CAD extractor exposes nested
    bounds.
    """

    try:
        from ..cad.cells import extract_cylindrical_solution, extract_explicit_cylindrical_solution

        cyl = extract_explicit_cylindrical_solution(expr, variables)
        if cyl is None:
            cyl = extract_cylindrical_solution(expr, variables, selected_only=True)
    except _EXPECTED_ERRORS:
        return None
    pieces: list[sp.Expr] = []
    n = len(tuple(variables))
    for cell in getattr(cyl, "cells", ()):
        if getattr(cell, "dimension", None) != n:
            continue
        levels = tuple(cell.levels)
        for i, level in enumerate(levels):
            for bound in (level.lower, level.upper):
                if bound in (-sp.oo, sp.oo):
                    continue
                constraints: list[sp.Expr] = []
                for j, other in enumerate(levels):
                    if i == j:
                        constraints.append(sp.Eq(other.variable, bound))
                    else:
                        constraints.append(other.as_formula(closed=True))
                pieces.append(sp.And(*constraints))
    if not pieces:
        return None
    return _simplify_region(sp.Or(*pieces))


def _prefer_vertical_boundary(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> bool:
    if len(variables) != 2:
        return False
    x_var, y_var = variables
    for atom in expr.atoms(sp.Equality):
        try:
            poly = sp.Poly(sp.expand(atom.lhs - atom.rhs), x_var, y_var)
        except _EXPECTED_ERRORS:
            continue
        if poly.degree(x_var) <= 1 and poly.degree(y_var) > 1:
            return True
    return False


def region_boundary(
    region: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> sp.Expr:
    """Return a formula for the boundary: closure(region) minus interior(region)."""
    expr = _normalize_formula(region)
    vars_ = _normalize_variables(variables, expr)
    if _syntactically_empty(expr):
        return sp.false
    if isinstance(expr, sp.Or):
        return _simplify_region(sp.Or(*[region_boundary(arg, vars_) for arg in expr.args]))
    if len(vars_) == 1:
        interval = _interval_from_piece(expr, vars_[0])
        if interval is not None:
            low, high, _, _ = interval
            endpoints = []
            if low is not -sp.oo:
                endpoints.append(sp.Eq(vars_[0], low))
            if high is not sp.oo and high != low:
                endpoints.append(sp.Eq(vars_[0], high))
            return _simplify_region(sp.Or(*endpoints) if endpoints else sp.false)
    explicit = _boundary_from_conjunction(expr, vars_)
    if explicit is not None and not _prefer_vertical_boundary(explicit, vars_):
        return explicit
    cad_vertical = _boundary_from_structured_cad_cells_2d(expr, vars_)
    if cad_vertical is not None:
        return cad_vertical
    vertical = _boundary_from_vertical_bounds_2d(expr, vars_)
    if vertical is not None:
        return vertical
    if explicit is not None:
        return explicit
    cad_nd = _boundary_from_cylindrical_solution_nd(expr, vars_)
    if cad_nd is not None:
        return cad_nd
    closure = region_closure(expr, vars_, strategy=strategy)
    interior = region_interior(expr, vars_, strategy=strategy)
    return _simplify_region(sp.And(closure, sp.Not(interior)))


def _atomic_equalities(expr: sp.Expr) -> list[sp.Expr]:
    if isinstance(expr, sp.And):
        out: list[sp.Expr] = []
        for arg in expr.args:
            out.extend(_atomic_equalities(arg))
        return out
    if isinstance(expr, sp.Equality):
        return [sp.simplify(expr.lhs - expr.rhs)]
    return []


def _dimension_from_equalities(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> int | None:
    eqs = _atomic_equalities(expr)
    if not eqs:
        return None
    linear_rows: list[list[sp.Expr]] = []
    nonlinear_count = 0
    for poly in eqs:
        try:
            p = sp.Poly(poly, *variables)
        except _EXPECTED_ERRORS:
            return None
        if p.total_degree() <= 1:
            linear_rows.append([p.coeff_monomial(var) for var in variables])
        else:
            nonlinear_count += 1
    if linear_rows and nonlinear_count == 0:
        return max(0, len(variables) - int(sp.Matrix(linear_rows).rank()))
    return max(0, len(variables) - nonlinear_count - (1 if linear_rows else 0))


def region_dimension(
    region: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> int:
    """Return a dimension estimate for common implicit semialgebraic regions."""
    expr = _normalize_formula(region)
    vars_ = _normalize_variables(variables, expr)
    if _syntactically_empty(expr):
        return -1
    interior = region_interior(expr, vars_, strategy=strategy)
    if not _syntactically_empty(interior):
        return len(vars_)
    estimate = _dimension_from_equalities(expr, vars_)
    if estimate is not None:
        return estimate
    if _is_relational(expr):
        return max(0, len(vars_) - 1)
    return 0


def _as_disjuncts(expr: sp.Expr) -> list[sp.Expr]:
    simplified = _simplify_region(expr)
    return list(simplified.args) if isinstance(simplified, sp.Or) else [simplified]


def _compare_endpoints(left: sp.Expr, right: sp.Expr) -> int:
    """Compare finite or infinite real interval endpoints exactly."""

    if left == right:
        return 0
    if left is -sp.oo or right is sp.oo:
        return -1
    if left is sp.oo or right is -sp.oo:
        return 1
    return compare_exact_reals(left, right)


def _interval_from_piece(piece: sp.Expr, variable: sp.Symbol):
    """Recover a one-dimensional interval description from a conjunction of bounds."""
    atoms = list(piece.args) if isinstance(piece, sp.And) else [piece]
    low = -sp.oo
    high = sp.oo
    low_closed = False
    high_closed = False
    for original in atoms:
        if not _is_relational(original):
            try:
                reduced = sp.reduce_inequalities([original], variable)
            except _EXPECTED_ERRORS:
                return None
            if reduced is sp.false or reduced == sp.false:
                return None
            if reduced != original:
                nested = _interval_from_piece(reduced, variable)
                if nested is None:
                    return None
                n_low, n_high, n_low_closed, n_high_closed = nested
                low_cmp = _compare_endpoints(n_low, low)
                if low_cmp > 0:
                    low, low_closed = n_low, n_low_closed
                elif low_cmp == 0:
                    low_closed = low_closed and n_low_closed
                high_cmp = _compare_endpoints(n_high, high)
                if high_cmp < 0:
                    high, high_closed = n_high, n_high_closed
                elif high_cmp == 0:
                    high_closed = high_closed and n_high_closed
                continue
            return None
        atom = original
        lhs, rhs = atom.lhs, atom.rhs
        if rhs == variable:
            lhs, rhs = rhs, lhs
            if isinstance(atom, sp.StrictGreaterThan):
                atom = sp.StrictLessThan(lhs, rhs)
            elif isinstance(atom, sp.GreaterThan):
                atom = sp.LessThan(lhs, rhs)
            elif isinstance(atom, sp.StrictLessThan):
                atom = sp.StrictGreaterThan(lhs, rhs)
            elif isinstance(atom, sp.LessThan):
                atom = sp.GreaterThan(lhs, rhs)
            else:
                atom = sp.Eq(lhs, rhs)
        if lhs != variable:
            return None
        if isinstance(atom, (sp.StrictGreaterThan, sp.GreaterThan)):
            cmp = _compare_endpoints(rhs, low)
            closed = isinstance(atom, sp.GreaterThan)
            if cmp > 0:
                low, low_closed = rhs, closed
            elif cmp == 0:
                low_closed = low_closed and closed
        elif isinstance(atom, (sp.StrictLessThan, sp.LessThan)):
            cmp = _compare_endpoints(rhs, high)
            closed = isinstance(atom, sp.LessThan)
            if cmp < 0:
                high, high_closed = rhs, closed
            elif cmp == 0:
                high_closed = high_closed and closed
        elif isinstance(atom, sp.Equality):
            low = high = rhs
            low_closed = high_closed = True
        else:
            return None
    cmp = _compare_endpoints(low, high)
    if cmp > 0 or (cmp == 0 and not (low_closed and high_closed)):
        return None
    return (low, high, low_closed, high_closed)


def _intervals_touch(left, right) -> bool:
    _, a_high, _, a_high_closed = left
    b_low, _, b_low_closed, _ = right
    cmp = _compare_endpoints(a_high, b_low)
    return cmp > 0 or (cmp == 0 and (a_high_closed or b_low_closed))


def _merge_intervals(left, right):
    """Merge overlapping/touching intervals without shrinking containment."""

    low, a_high, low_closed, a_high_closed = left
    _, b_high, _, b_high_closed = right
    cmp = _compare_endpoints(a_high, b_high)
    if cmp > 0:
        return left
    if cmp < 0:
        return (low, b_high, low_closed, b_high_closed)
    return (low, a_high, low_closed, a_high_closed or b_high_closed)


def _piece_from_interval(interval, variable: sp.Symbol) -> sp.Expr:
    low, high, low_closed, high_closed = interval
    clauses: list[sp.Expr] = []
    if low is not -sp.oo:
        clauses.append(variable >= low if low_closed else variable > low)
    if high is not sp.oo:
        clauses.append(variable <= high if high_closed else variable < high)
    return _simplify_region(sp.And(*clauses) if clauses else sp.true)


def _components_1d(expr: sp.Expr, variable: sp.Symbol) -> tuple[sp.Expr, ...] | None:
    intervals = []
    for piece in _as_disjuncts(expr):
        interval = _interval_from_piece(piece, variable)
        if interval is None:
            return None
        intervals.append(interval)
    intervals.sort(key=cmp_to_key(lambda a, b: _compare_endpoints(a[0], b[0])))
    merged = []
    for interval in intervals:
        if not merged:
            merged.append(interval)
            continue
        prev = merged[-1]
        if _intervals_touch(prev, interval):
            merged[-1] = _merge_intervals(prev, interval)
        else:
            merged.append(interval)
    return tuple(_piece_from_interval(interval, variable) for interval in merged)


def region_components(
    region: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> tuple[sp.Expr, ...]:
    """Return connected-component formulas for simple explicit cases.

    The first implementation is exact for one-dimensional semialgebraic sets
    reducible by SymPy's inequality reducer. For higher-dimensional explicit
    disjunctions it returns the top-level nonempty pieces without claiming a
    full CAD adjacency computation.
    """
    expr = _normalize_formula(region)
    vars_ = _normalize_variables(variables, expr)
    if _syntactically_empty(expr):
        return ()
    if len(vars_) == 1:
        components = _components_1d(expr, vars_[0])
        if components is not None:
            return components
    return tuple(piece for piece in _as_disjuncts(expr) if not _syntactically_empty(piece))
