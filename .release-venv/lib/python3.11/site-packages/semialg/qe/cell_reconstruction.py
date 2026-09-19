"""Exact formula reconstruction from projected CAD cells."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..cad_algorithms.lifting.stack import CADCell
from ..reconstruct.merge import compressed_formula_from_cells, dnf_formula_from_cells
from ..simplify.boolean import simplify_boolean


def finite_variety_formula(cells: Sequence[CADCell], variables: Sequence[sp.Symbol]) -> sp.Expr:
    """Reconstruct a finite variety sub-CAD as a union of exact points."""

    if not cells:
        return sp.false
    terms: list[sp.Expr] = []
    for cell in cells:
        if len(cell.sample) < len(variables):
            continue
        terms.append(
            sp.And(
                *(
                    sp.Eq(variable, sample_to_expr(cell.sample[index]))
                    for index, variable in enumerate(variables)
                )
            )
        )
    return sp.Or(*terms) if terms else sp.false


def cells_to_formula(
    cells: Sequence[CADCell],
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    form: str = "nested",
    max_terms: int = 512,
) -> sp.Expr:
    """Convert a cell union to a conservative exact SymPy formula.

    ``form="nested"`` shares common cylindrical prefixes and merges contiguous
    stack blocks. ``form="dnf"`` emits one disjunct per selected cell.
    """

    if not cells:
        return sp.false
    if form == "dnf":
        return simplify_boolean(dnf_formula_from_cells(cells, variables, cells_by_level))
    result = compressed_formula_from_cells(cells, variables, cells_by_level, max_terms=max_terms)
    if result.stats.fallback_used:
        return simplify_boolean(result.formula)
    return result.formula


__all__ = ["cells_to_formula", "finite_variety_formula"]
