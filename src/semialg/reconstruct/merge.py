from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from ..cad.lifting.stack import CADCell
from .cylindrical import path_condition
from .nested import NestedFormulaResult, nested_formula_from_cells


def dnf_formula_from_cells(
    cells: Sequence[CADCell],
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    closed: bool = False,
) -> sp.Expr:
    """Conservative one-disjunct-per-cell fallback formula."""

    if not cells:
        return sp.false
    return sp.Or(
        *(path_condition(cell, variables, cells_by_level, closed=closed) for cell in cells)
    )


def compressed_formula_from_cells(
    cells: Sequence[CADCell],
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    max_terms: int = 512,
    closed: bool = False,
) -> NestedFormulaResult:
    """Return a nested formula, falling back to DNF if compression overflows."""

    result = nested_formula_from_cells(
        cells, variables, cells_by_level, max_terms=max_terms, closed=closed
    )
    if not result.stats.fallback_used:
        return result
    fallback = dnf_formula_from_cells(cells, variables, cells_by_level, closed=closed)
    return NestedFormulaResult(fallback, result.stats)


__all__ = ["compressed_formula_from_cells", "dnf_formula_from_cells"]
