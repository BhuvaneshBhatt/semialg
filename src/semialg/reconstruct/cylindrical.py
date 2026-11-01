from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from ..algebraic.comparison import compare_samples
from ..algebraic.roots import isolate_real_roots
from ..algebraic.samples import AlgebraicRoot, sample_to_expr
from ..cad_algorithms.bounds import (
    AlgebraicRootFunction,
    CADBound,
    DelineabilityCertificate,
    as_cad_bound,
)
from ..cad_algorithms.lifting.stack import CADCell
from ..cad_algorithms.polynomial_utils import exact_univariate_poly
from ..errors import ExactEvaluationFailure
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


def _section_root_index(cell: CADCell, variable: sp.Symbol, base_variables) -> int:
    """Identify a section among the distinct roots of its defining polynomial.

    A stack can interleave roots from several polynomials. Its section position
    is therefore not a family-local root rank. Reuse the lifting root's rank
    when its provenance matches; otherwise isolate the specialized family and
    identify the sample by exact comparison, including algebraic coefficients.
    """
    polynomial = cell.section_polynomial
    sample = cell.sample[cell.level - 1]
    if isinstance(sample, AlgebraicRoot):
        defining = (
            sample.fiber_context.polynomial
            if sample.fiber_context is not None
            else sample.polynomial
        )
        if sp.expand(polynomial - defining.as_expr()) == 0:
            return sample.root_index
    substitutions = {var: sample_to_expr(cell.sample[i]) for i, var in enumerate(base_variables)}
    specialized = sp.expand(polynomial.subs(substitutions))
    if (
        isinstance(sample, AlgebraicRoot)
        and sp.expand(specialized - sample.polynomial.as_expr()) == 0
    ):
        return sample.root_index
    roots = isolate_real_roots(
        exact_univariate_poly(specialized, variable, algebraic_extension=False)
    )
    for index, root in enumerate(roots):
        if compare_samples(root, sample) == 0:
            return index
    raise ExactEvaluationFailure("could not identify the section in its defining real-root family")


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
            certificate.root_index
            if certificate is not None
            else _section_root_index(cell, variable, base_variables)
        )
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
