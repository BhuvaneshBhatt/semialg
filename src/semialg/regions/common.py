from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from ..model import Cell, QEResult


def cell_truth_list(obj):
    if isinstance(obj, QEResult):
        return qe_cells_and_vars(obj)
    return list(obj), None


def cell_to_formula(
    cell: Cell,
    variables: Sequence[sp.Symbol],
    *,
    left_closed: bool,
    right_closed: bool,
) -> sp.Expr:
    clauses = []
    for index, variable in enumerate(variables[: cell.level]):
        left, right = cell.intervals[index]
        point_cell = left is not None and right is not None and sp.simplify(left - right) == 0
        if point_cell:
            if left_closed and right_closed:
                clauses.append(sp.Eq(variable, left))
            else:
                return sp.false
            continue
        if left is not None:
            clauses.append(variable >= left if left_closed else variable > left)
        if right is not None:
            clauses.append(variable <= right if right_closed else variable < right)
    return sp.And(*clauses) if clauses else sp.true


def region_formula(
    cells_with_truth, variables, *, left_closed: bool, right_closed: bool
) -> sp.Expr:
    pieces = []
    for cell, truth in cells_with_truth:
        if not truth:
            continue
        piece = cell_to_formula(cell, variables, left_closed=left_closed, right_closed=right_closed)
        if piece is not sp.false:
            pieces.append(piece)
    if not pieces:
        return sp.false
    return sp.simplify_logic(sp.Or(*pieces), form="dnf")


def qe_cells_and_vars(qe_result: QEResult):
    if qe_result.free_level_cells is None:
        raise ValueError("Region operations require free-level cells, not a sentence result")
    return list(qe_result.free_level_cells), qe_result.free_vars
