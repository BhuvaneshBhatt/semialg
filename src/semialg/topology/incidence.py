from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..cad.decomposition import CompleteCAD
from ..cad.lifting.stack import CADCell
from ..reconstruct.cylindrical import path_condition


def cell_dimension(cell: CADCell, cells_by_level: Mapping[int, Sequence[CADCell]]) -> int:
    """Return the Euclidean dimension of a final CAD cell.

    A sector contributes one dimension and a section contributes zero. The
    dimension is read from the complete chain of ancestors, not merely from the
    final fiber cell, because a two-dimensional section over a one-dimensional
    base is a curve, while a section over a point is a point.
    """

    dim = 0
    for level in range(1, cell.level + 1):
        ancestor = cell_at_index(cells_by_level, cell.index[:level])
        if ancestor.kind == "sector":
            dim += 1
    return dim


def cell_at_index(
    cells_by_level: Mapping[int, Sequence[CADCell]], index: tuple[int, ...]
) -> CADCell:
    level = len(index)
    for cell in cells_by_level[level]:
        if cell.index == index:
            return cell
    raise KeyError(index)


def final_cells(cad: CompleteCAD) -> tuple[CADCell, ...]:
    return cad.cells_by_level.get(len(cad.tower.variables), tuple())


def cell_sample_subs(cell: CADCell, variables: Sequence[sp.Symbol]) -> dict[sp.Symbol, sp.Expr]:
    return {var: sample_to_expr(sample) for var, sample in zip(variables, cell.sample, strict=True)}


def is_cell_in_closure(
    target: CADCell,
    source: CADCell,
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
) -> bool:
    """Return whether ``target`` is contained in the closure of ``source``.

    CAD closure/incidence is recursive and can become subtle when fiber root
    functions degenerate over boundary base cells. The implementation uses the
    cylindrical closed path formula for ``source`` and tests it on the sample of
    ``target``. Because the CAD is sign-invariant for all projection/root
    boundary polynomials, this sample test is exact for the cell-level closure
    operations used by this package. A dimension guard prevents sectors from
    being treated as contained in the closure of lower-dimensional sections.
    """

    if target.level != source.level:
        return False
    if cell_dimension(target, cells_by_level) > cell_dimension(source, cells_by_level):
        return False
    if target.index == source.index:
        return True
    condition = path_condition(source, variables, cells_by_level, closed=True)
    return _truth_at_cell_sample(condition, target, variables)


def _truth_at_cell_sample(
    condition: sp.Expr, cell: CADCell, variables: Sequence[sp.Symbol]
) -> bool:
    value = sp.simplify(condition.subs(cell_sample_subs(cell, variables)))
    if value is sp.true or value == sp.true:
        return True
    if value is sp.false or value == sp.false:
        return False
    try:
        return bool(value)
    except TypeError:
        pass
    try:
        return bool(sp.N(value))
    except Exception:
        return False


def closures_intersect(
    left: CADCell,
    right: CADCell,
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    candidates: Sequence[CADCell] | None = None,
) -> bool:
    """Return whether two CAD cell closures intersect in the final CAD level."""

    if left.index == right.index:
        return True
    ambient = candidates or tuple(cells_by_level.get(left.level, ()))
    max_dim = min(cell_dimension(left, cells_by_level), cell_dimension(right, cells_by_level))
    for cell in ambient:
        if cell_dimension(cell, cells_by_level) > max_dim:
            continue
        if is_cell_in_closure(cell, left, variables, cells_by_level) and is_cell_in_closure(
            cell, right, variables, cells_by_level
        ):
            return True
    return False


__all__ = [
    "cell_at_index",
    "cell_dimension",
    "cell_sample_subs",
    "closures_intersect",
    "final_cells",
    "is_cell_in_closure",
]
