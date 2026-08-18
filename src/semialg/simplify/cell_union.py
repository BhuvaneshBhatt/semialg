from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..reconstruct.cylindrical import path_condition
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
    return sp.simplify_logic(sp.Or(*pieces), form="dnf")


def cell_union_to_formula(cell_union: CellUnion) -> sp.Expr:
    if len(cell_union.variables) == 0:
        return sp.true if cell_union.cells else sp.false
    if len(cell_union.variables) == 1:
        return intervals_to_formula(cell_union.variables[0], cell_union_to_intervals(cell_union))
    return _multivariate_cell_form(cell_union)


__all__ = ["cell_to_interval_1d", "cell_union_to_intervals", "cell_union_to_formula"]
