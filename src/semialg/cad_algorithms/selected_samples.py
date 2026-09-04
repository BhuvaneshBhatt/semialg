"""Lightweight exact sample extraction from selected terminal CAD cells."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..formula import formula_polynomials, parse_formula
from ..normalization import normalize_formula, normalize_problem_variables
from .decomposition import decomp_collins_complete


@dataclass(frozen=True)
class SelectedCADSample:
    index: tuple[int, ...]
    point: Mapping[sp.Symbol, sp.Expr]


def extract_selected_cad_samples(
    condition: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    limit: int | None = None,
) -> tuple[SelectedCADSample, ...]:
    """Return exact representative points without reconstructing cylindrical cells."""
    formula = normalize_formula(condition)
    vars_ = normalize_problem_variables(variables, formula)
    parsed = parse_formula(formula)
    from ..qe.complete import evaluate_formula_on_cell

    cad_obj = decomp_collins_complete(formula_polynomials(parsed), vars_)
    result: list[SelectedCADSample] = []
    for cell in cad_obj.cells:
        if evaluate_formula_on_cell(parsed, cell, vars_):
            point = {
                var: sample_to_expr(value) for var, value in zip(vars_, cell.sample, strict=True)
            }
            result.append(SelectedCADSample(cell.index, point))
            if limit is not None and len(result) >= limit:
                break
    return tuple(result)


__all__ = ["SelectedCADSample", "extract_selected_cad_samples"]
