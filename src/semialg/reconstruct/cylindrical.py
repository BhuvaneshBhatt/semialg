from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from .._zero_testing import certified_equal
from ..algebraic.samples import sample_to_expr
from ..cad_algorithms.bounds import (
    AlgebraicRootFunction,
    CADBound,
    DelineabilityCertificate,
    as_cad_bound,
)
from ..cad_algorithms.lifting.stack import CADCell
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


def section_value_bound(
    cell: CADCell,
    variable: sp.Symbol,
    *,
    base_variables: Sequence[sp.Symbol] = (),
    certificate: DelineabilityCertificate | None = None,
    closed: bool = True,
) -> CADBound:
    """Return a typed value for a section cell.

    Variable-dependent algebraic sections are represented by
    :class:`AlgebraicRootFunction`; ordinary fixed algebraic/rational sections
    retain their exact sample representation.
    """

    poly = cell.section_polynomial
    if (
        poly is not None
        and cell.root_index is not None
        and variable in sp.sympify(poly).free_symbols
    ):
        coefficient_symbols = sp.sympify(poly).free_symbols - {variable}
        local_root_index = int(
            certificate.root_index if certificate is not None else cell.root_index
        )
        if certificate is None:
            try:
                base_subs = {
                    base_variables[i]: sample_to_expr(cell.sample[i])
                    for i in range(min(cell.level - 1, len(base_variables)))
                }
                specialized = sp.expand(poly.subs(base_subs))
                exact_poly = sp.Poly(specialized, variable, extension=True)
                roots = tuple(root for root in exact_poly.all_roots() if root.is_real is not False)
                sample_value = sample_to_expr(cell.sample[cell.level - 1])
                matches = [
                    i for i, root in enumerate(roots) if certified_equal(root, sample_value) is True
                ]
                if len(matches) == 1:
                    local_root_index = matches[0]
            except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
                pass
        if coefficient_symbols:
            return AlgebraicRootFunction(
                polynomial=sp.expand(poly),
                fiber_variable=variable,
                root_index=local_root_index,
                base_variables=tuple(base_variables),
                base_index=cell.parent_index,
                certificate=certificate,
                stack_root_index=cell.root_index,
                closed=closed,
            )
        expr = fiber_root_expr(poly, variable, local_root_index)
        if expr is not None:
            return as_cad_bound(expr, closed=closed)
    point = cell.lower_bound or cell.upper_bound
    if point is None:
        raise ValueError("section cell has no reconstructible value")
    return as_cad_bound(point, closed=closed)


def section_value_expr(
    cell: CADCell, variable: sp.Symbol, *, base_variables: Sequence[sp.Symbol] = ()
) -> sp.Expr:
    """Return the symbolic value of a section cell as a radical/root function."""

    return section_value_bound(cell, variable, base_variables=base_variables).as_expr()


def _sector_bound_expr(
    sector: CADCell,
    variable: sp.Symbol,
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    side: str,
    base_variables: Sequence[sp.Symbol] = (),
) -> sp.Expr | None:
    offset = -1 if side == "left" else 1
    section = _section_at_position(
        cells_by_level, sector.level, sector.parent_index, sector.stack_position + offset
    )
    if section is not None:
        return section_value_expr(section, variable, base_variables=base_variables)
    bound = sector.lower_bound if side == "left" else sector.upper_bound
    return None if bound is None else sample_to_expr(bound)


def level_cell_condition(
    cell: CADCell,
    variable: sp.Symbol,
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    closed: bool = False,
    base_variables: Sequence[sp.Symbol] = (),
) -> sp.Expr:
    """Return a symbolic condition for one level of a CAD path."""

    if cell.kind == "section":
        return sp.Eq(variable, section_value_expr(cell, variable, base_variables=base_variables))
    left = _sector_bound_expr(
        cell, variable, cells_by_level, side="left", base_variables=base_variables
    )
    right = _sector_bound_expr(
        cell, variable, cells_by_level, side="right", base_variables=base_variables
    )
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
            level_cell_condition(
                level_cell,
                variables[level - 1],
                cells_by_level,
                closed=closed,
                base_variables=variables[: level - 1],
            )
        )
    kept = [piece for piece in pieces if piece is not sp.true and piece != sp.true]
    return sp.And(*kept) if kept else sp.true


__all__ = ["level_cell_condition", "path_condition", "section_value_bound", "section_value_expr"]
