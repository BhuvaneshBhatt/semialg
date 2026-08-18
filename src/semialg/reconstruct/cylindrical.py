from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..cad.lifting.stack import CADCell
from .radicals import fiber_root_expr


def _cells_for_parent(
    cells_by_level: Mapping[int, Sequence[CADCell]], level: int, parent: tuple[int, ...] | None
) -> tuple[CADCell, ...]:
    return tuple(
        sorted(
            (cell for cell in cells_by_level.get(level, ()) if cell.parent_index == parent),
            key=lambda cell: cell.stack_position,
        )
    )


def _section_at_position(
    cells_by_level: Mapping[int, Sequence[CADCell]],
    level: int,
    parent: tuple[int, ...] | None,
    position: int,
) -> CADCell | None:
    for cell in _cells_for_parent(cells_by_level, level, parent):
        if cell.stack_position == position and cell.kind == "section":
            return cell
    return None


def section_value_expr(cell: CADCell, variable: sp.Symbol) -> sp.Expr:
    """Return the symbolic value of a section cell as a radical/root function."""

    expr = fiber_root_expr(cell.section_polynomial, variable, cell.root_index)
    if expr is not None:
        return expr
    point = cell.lower_bound or cell.upper_bound
    return sample_to_expr(point)


def _sector_bound_expr(
    sector: CADCell,
    variable: sp.Symbol,
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    side: str,
) -> sp.Expr | None:
    offset = -1 if side == "left" else 1
    section = _section_at_position(
        cells_by_level, sector.level, sector.parent_index, sector.stack_position + offset
    )
    if section is not None:
        return section_value_expr(section, variable)
    bound = sector.lower_bound if side == "left" else sector.upper_bound
    return None if bound is None else sample_to_expr(bound)


def level_cell_condition(
    cell: CADCell,
    variable: sp.Symbol,
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    closed: bool = False,
) -> sp.Expr:
    """Return a symbolic condition for one level of a CAD path."""

    if cell.kind == "section":
        return sp.Eq(variable, section_value_expr(cell, variable))
    left = _sector_bound_expr(cell, variable, cells_by_level, side="left")
    right = _sector_bound_expr(cell, variable, cells_by_level, side="right")
    pieces: list[sp.Expr] = []
    if left is not None:
        pieces.append(variable >= left if closed else variable > left)
    if right is not None:
        pieces.append(variable <= right if closed else variable < right)
    return sp.And(*pieces) if pieces else sp.true


def path_condition(
    cell: CADCell,
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    closed: bool = False,
) -> sp.Expr:
    """Return a cylindrical formula for the path ending at ``cell``."""

    pieces: list[sp.Expr] = []
    for level in range(1, cell.level + 1):
        prefix = cell.index[:level]
        level_cell = next(item for item in cells_by_level[level] if item.index == prefix)
        pieces.append(
            level_cell_condition(level_cell, variables[level - 1], cells_by_level, closed=closed)
        )
    kept = [piece for piece in pieces if piece is not sp.true and piece != sp.true]
    return sp.And(*kept) if kept else sp.true


__all__ = ["level_cell_condition", "path_condition", "section_value_expr"]
