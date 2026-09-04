from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..cad_algorithms.polynomial_utils import polynomial_key
from ..reconstruct.cylindrical import path_condition
from .boolean import simplify_boolean
from .intervals import Interval1D, intervals_to_formula, merge_intervals

if TYPE_CHECKING:  # pragma: no cover
    from ..qe.complete import CellUnion


def cell_to_interval_1d(cell) -> Interval1D:
    """Convert a one-dimensional CAD cell to an interval."""

    left, right = cell.interval or (None, None)
    if cell.kind == "section":
        point = sample_to_expr(left if left is not None else right)
        return Interval1D(point, point, True, True)
    return Interval1D(
        sample_to_expr(left) if left is not None else None,
        sample_to_expr(right) if right is not None else None,
        False,
        False,
    )


def cell_union_to_intervals(cell_union: CellUnion) -> tuple[Interval1D, ...]:
    if len(cell_union.variables) != 1:
        raise ValueError("interval reconstruction requires exactly one free variable")
    return merge_intervals(tuple(cell_to_interval_1d(cell) for cell in cell_union.cells))


def _children(
    cells_by_level: Mapping[int, Sequence[object]], parent_index: tuple[int, ...], child_level: int
) -> tuple[object, ...]:
    return tuple(
        cell for cell in cells_by_level.get(child_level, ()) if cell.index[:-1] == parent_index
    )


def compress_sel_idxs(
    selected: set[tuple[int, ...]],
    *,
    free_level: int,
    cells_by_level: Mapping[int, Sequence[object]],
) -> set[tuple[int, ...]]:
    """Compress selected final cells up the cylindrical tree when safe.

    If every child in a stack over a parent is selected, the whole stack is
    replaced by the parent cell. Repeating from top to bottom gives a compact
    cylindrical cover. This is conservative: it never invents non-cylindrical
    implications and never merges partial stacks.
    """
    current = set(selected)
    for level in range(free_level, 1, -1):
        parent_to_selected: dict[tuple[int, ...], set[tuple[int, ...]]] = {}
        for index in tuple(current):
            if len(index) == level:
                parent_to_selected.setdefault(index[:-1], set()).add(index)
        for parent, chosen_children in parent_to_selected.items():
            all_children = {child.index for child in _children(cells_by_level, parent, level)}
            if all_children and chosen_children == all_children:
                current.difference_update(chosen_children)
                current.add(parent)
    return current


def _multivariate_cell_form(cell_union: CellUnion) -> sp.Expr:
    cells_by_level = getattr(cell_union, "cells_by_level", None) or {}
    free_level = len(cell_union.variables)
    if not cells_by_level or free_level <= 1:
        return cell_union.formula
    selected = {cell.index for cell in cell_union.cells}
    compressed = compress_sel_idxs(selected, free_level=free_level, cells_by_level=cells_by_level)
    pieces: list[sp.Expr] = []
    for index in sorted(compressed):
        level = len(index)
        try:
            cell = next(item for item in cells_by_level[level] if item.index == index)
        except StopIteration:
            continue
        pieces.append(path_condition(cell, cell_union.variables, cells_by_level))
    if not pieces:
        return sp.false
    return simplify_boolean(sp.Or(*pieces))


def _sign_relation(poly: sp.Expr, signs: frozenset[int]) -> sp.Expr | None:
    """Return the canonical polynomial relation for a selected sign subset."""

    if signs == frozenset({-1}):
        return poly < 0
    if signs == frozenset({0}):
        return sp.Eq(poly, 0)
    if signs == frozenset({1}):
        return poly > 0
    if signs == frozenset({-1, 0}):
        return poly <= 0
    if signs == frozenset({0, 1}):
        return poly >= 0
    if signs == frozenset({-1, 1}):
        return sp.Ne(poly, 0)
    if signs == frozenset({-1, 0, 1}):
        return sp.true
    if not signs:
        return sp.false
    return None


def polynomial_relation_from_cell_union(cell_union: CellUnion) -> sp.Expr | None:
    """Recover a single polynomial sign relation when it exactly selects the cells.

    This is an exact reconstruction shortcut, not a heuristic.  Every
    full-dimensional free-variable CAD cell has a certified sign for each
    projection polynomial active on its ancestor level.  A candidate relation
    is returned only when membership in the selected cell union is completely
    determined by one polynomial's sign over *all* free-variable cells.

    The result is therefore directly reusable by polynomial CAD/QE and avoids
    opaque root-function/``Abs`` boundary formulas when a simpler sign
    description exists.
    """

    free_level = len(cell_union.variables)
    cells_by_level = cell_union.cells_by_level
    projection_polynomials = getattr(cell_union, "projection_polynomials", None) or {}
    all_cells = tuple(cells_by_level.get(free_level, ()))
    if free_level == 0 or not all_cells or not projection_polynomials:
        return None

    selected = {cell.index for cell in cell_union.cells}
    if not selected:
        return sp.false
    if selected == {cell.index for cell in all_cells}:
        return sp.true

    candidates: list[sp.Expr] = []
    for level in range(free_level, 0, -1):
        for expr in projection_polynomials.get(level, ()):
            key = polynomial_key(sp.Poly(expr, *cell_union.variables[:level], domain="EX"))
            by_sign: dict[int, set[bool]] = {}
            usable = True
            for leaf in all_cells:
                ancestor_index = leaf.index[:level]
                try:
                    ancestor = next(
                        cell
                        for cell in cells_by_level.get(level, ())
                        if cell.index == ancestor_index
                    )
                except StopIteration:
                    usable = False
                    break
                sign = ancestor.signs.get(key)
                if sign not in (-1, 0, 1):
                    usable = False
                    break
                by_sign.setdefault(int(sign), set()).add(leaf.index in selected)
            if not usable or any(len(values) != 1 for values in by_sign.values()):
                continue
            selected_signs = frozenset(sign for sign, values in by_sign.items() if True in values)
            relation = _sign_relation(sp.expand(expr), selected_signs)
            if relation is not None:
                candidates.append(relation)

    if not candidates:
        return None
    return min(
        candidates,
        key=lambda expr: (sp.count_ops(expr), len(sp.sstr(expr)), sp.default_sort_key(expr)),
    )


def cell_union_to_formula(cell_union: CellUnion) -> sp.Expr:
    polynomial = polynomial_relation_from_cell_union(cell_union)
    if polynomial is not None:
        return polynomial
    if len(cell_union.variables) == 0:
        return sp.true if cell_union.cells else sp.false
    if len(cell_union.variables) == 1:
        return intervals_to_formula(cell_union.variables[0], cell_union_to_intervals(cell_union))
    return _multivariate_cell_form(cell_union)


__all__ = [
    "cell_to_interval_1d",
    "cell_union_to_intervals",
    "cell_union_to_formula",
    "polynomial_relation_from_cell_union",
]
