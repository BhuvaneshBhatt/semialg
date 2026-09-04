from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..cad_algorithms.decomposition import CompleteCAD
from ..cad_algorithms.lifting.stack import CADCell
from ..reconstruct.merge import compressed_formula_from_cells
from ..simplify.boolean import simplify_boolean
from ..simplify.formula import simplify_qe_formula
from .incidence import cell_dimension, final_cells, is_cell_in_closure


@dataclass(frozen=True)
class TopologyResult:
    operation: str
    cells: tuple[CADCell, ...]
    formula: sp.Expr


def closure_cells(
    selected: Iterable[CADCell],
    cad: CompleteCAD,
    variables: Sequence[sp.Symbol],
) -> tuple[CADCell, ...]:
    """Return final CAD cells contained in the Euclidean closure."""

    selected_tuple = tuple(selected)
    ambient = final_cells(cad)
    cells_by_level = cad.cells_by_level
    if len(tuple(variables)) == 1:
        selected_pos = {cell.stack_position for cell in selected_tuple}
        closed = [
            cell
            for cell in ambient
            if cell.stack_position in selected_pos
            or (
                cell.kind == "section"
                and (
                    cell.stack_position - 1 in selected_pos
                    or cell.stack_position + 1 in selected_pos
                )
            )
        ]
        return tuple(sorted(closed, key=lambda cell: cell.index))
    closed = [
        cell
        for cell in ambient
        if any(
            is_cell_in_closure(cell, source, variables, cells_by_level) for source in selected_tuple
        )
    ]
    return tuple(sorted(closed, key=lambda cell: cell.index))


def interior_cells(
    selected: Iterable[CADCell],
    cad: CompleteCAD,
    variables: Sequence[sp.Symbol],
) -> tuple[CADCell, ...]:
    """Return final CAD cells contained in the Euclidean interior.

    Interior is computed semantically as ``ambient - closure(ambient - S)``.
    This is deliberately *not* the same as keeping only selected
    full-dimensional CAD cells: an adapted CAD may split an open region by a
    lower-dimensional section that is itself contained in the interior.  For
    example, ``[0, 1] | [1, 2]`` contains the section ``x = 1`` as an interior
    point of the union even though that section is zero-dimensional.
    """

    selected_tuple = tuple(selected)
    selected_indices = {cell.index for cell in selected_tuple}
    ambient = final_cells(cad)
    complement = tuple(cell for cell in ambient if cell.index not in selected_indices)
    complement_closure = {cell.index for cell in closure_cells(complement, cad, variables)}
    return tuple(
        sorted(
            (cell for cell in selected_tuple if cell.index not in complement_closure),
            key=lambda cell: cell.index,
        )
    )


def boundary_cells(
    selected: Iterable[CADCell],
    cad: CompleteCAD,
    variables: Sequence[sp.Symbol],
) -> tuple[CADCell, ...]:
    """Return the Euclidean boundary using CAD closure semantics.

    The boundary is ``closure(S) ∩ closure(complement(S))``.  Expressing it
    this way avoids syntactic boundary artifacts at internal CAD sections.
    """

    selected_tuple = tuple(selected)
    selected_indices = {cell.index for cell in selected_tuple}
    ambient = final_cells(cad)
    complement = tuple(cell for cell in ambient if cell.index not in selected_indices)
    closed_selected = {cell.index for cell in closure_cells(selected_tuple, cad, variables)}
    closed_complement = {cell.index for cell in closure_cells(complement, cad, variables)}
    boundary_indices = closed_selected & closed_complement
    return tuple(
        sorted(
            (cell for cell in ambient if cell.index in boundary_indices),
            key=lambda cell: cell.index,
        )
    )


def exterior_cells(
    selected: Iterable[CADCell],
    cad: CompleteCAD,
    variables: Sequence[sp.Symbol],
) -> tuple[CADCell, ...]:
    """Return the Euclidean exterior: interior of the complement."""

    closure = {cell.index for cell in closure_cells(selected, cad, variables)}
    dim = len(tuple(variables))
    return tuple(
        sorted(
            (
                cell
                for cell in final_cells(cad)
                if cell.index not in closure and cell_dimension(cell, cad.cells_by_level) == dim
            ),
            key=lambda cell: cell.index,
        )
    )


def cells_formula(
    cells: Iterable[CADCell],
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    closed: bool = False,
    max_terms: int = 512,
) -> sp.Expr:
    selected = tuple(cells)
    if not selected:
        return sp.false
    result = compressed_formula_from_cells(
        selected, variables, cells_by_level, closed=closed, max_terms=max_terms
    )
    formula = result.formula
    if result.stats.fallback_used:
        formula = simplify_boolean(formula)
    return simplify_qe_formula(formula, implication_minimize=False)


def apply_topological_operation(
    selected: Iterable[CADCell],
    cad: CompleteCAD,
    variables: Sequence[sp.Symbol],
    operation: str | None,
) -> TopologyResult:
    if operation is None:
        cells = tuple(sorted(selected, key=lambda cell: cell.index))
    elif operation == "closure":
        cells = closure_cells(selected, cad, variables)
    elif operation == "interior":
        cells = interior_cells(selected, cad, variables)
    elif operation == "boundary":
        cells = boundary_cells(selected, cad, variables)
    elif operation == "exterior":
        cells = exterior_cells(selected, cad, variables)
    else:
        raise ValueError(f"unsupported topological operation: {operation!r}")
    return TopologyResult(
        operation or "identity", cells, cells_formula(cells, variables, cad.cells_by_level)
    )


__all__ = [
    "TopologyResult",
    "apply_topological_operation",
    "boundary_cells",
    "cells_formula",
    "closure_cells",
    "exterior_cells",
    "interior_cells",
]
