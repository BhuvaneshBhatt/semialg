from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ..cad.decomposition import CompleteCAD
from ..cad.lifting.stack import CADCell
from ..reconstruct.cylindrical import path_condition
from ..simplify.result import simplify_qe_formula
from ..topology.operations import boundary_cells, cells_formula, interior_cells
from .exceptional import (
    ExceptionalCause,
    exceptional_formula_from_causes,
    relevant_causes_for_cells,
)


@dataclass(frozen=True)
class GenericSplit:
    """Cell-level split used by generic CAD output."""

    generic_cells: tuple[CADCell, ...]
    exceptional_cells: tuple[CADCell, ...]
    exceptional_polys: tuple[sp.Poly, ...]
    exceptional_causes: tuple[ExceptionalCause, ...]
    generic_formula: sp.Expr
    exceptional_formula: sp.Expr


def _generic_cell_formula(
    cells: Sequence[CADCell],
    variables: Sequence[sp.Symbol],
    cells_by_level,
) -> sp.Expr:
    pieces = [path_condition(cell, variables, cells_by_level, closed=False) for cell in cells]
    kept = [piece for piece in pieces if piece is not sp.false and piece != sp.false]
    if not kept:
        return sp.false
    return simplify_qe_formula(sp.Or(*kept), implication_minimize=False)


def generic_split_from_cells(
    selected_cells: Sequence[CADCell],
    cad: CompleteCAD,
    variables: Sequence[sp.Symbol],
    causes: Sequence[ExceptionalCause] = (),
) -> GenericSplit:
    """Split a selected CAD cell set into full-dimensional and exceptional parts."""

    vars_tuple = tuple(variables)
    generic_cells = interior_cells(selected_cells, cad, vars_tuple)
    exceptional_cells = boundary_cells(selected_cells, cad, vars_tuple)
    generic_formula = _generic_cell_formula(generic_cells, vars_tuple, cad.cells_by_level)
    cell_exceptional = cells_formula(exceptional_cells, vars_tuple, cad.cells_by_level)
    relevant = relevant_causes_for_cells(causes, cell_exceptional)
    cause_exceptional = exceptional_formula_from_causes(relevant)
    exceptional_formula = cause_exceptional if cause_exceptional != sp.false else cell_exceptional
    return GenericSplit(
        generic_cells=tuple(generic_cells),
        exceptional_cells=tuple(exceptional_cells),
        exceptional_polys=tuple(cause.polynomial for cause in relevant),
        exceptional_causes=tuple(relevant),
        generic_formula=simplify_qe_formula(generic_formula, implication_minimize=False),
        exceptional_formula=simplify_qe_formula(exceptional_formula, implication_minimize=False),
    )


__all__ = ["GenericSplit", "generic_split_from_cells"]
