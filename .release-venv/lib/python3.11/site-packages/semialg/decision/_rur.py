from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from .._errors import EXACT_OPERATION_ERRORS as _RECOVERABLE_RUR_ERRORS
from ..algebraic.rational_univariate import solve_formula_with_rur


def finite_points_formula(
    points: Sequence[Mapping[sp.Symbol, sp.Expr]],
    variables: Sequence[sp.Symbol],
) -> sp.Expr:
    """Return an exact finite-set formula from point assignments."""

    if not points:
        return sp.false
    formulas = [
        sp.And(*(sp.Eq(var, sp.sympify(point[var])) for var in variables)) for point in points
    ]
    return sp.Or(*formulas) if len(formulas) > 1 else formulas[0]


def try_rur_formula(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    max_solutions: int | None = None,
):
    """Try the exact RUR backend for a supported finite Boolean formula.

    ``None`` means that the formula is outside the supported fragment or that
    the backend did not certify a complete answer.  In particular, an unknown
    or assignment-free partial result is never interpreted as unsatisfiable.
    """

    if not variables:
        return None
    try:
        result = solve_formula_with_rur(
            formula,
            tuple(variables),
            real=True,
            max_solutions=max_solutions,
        )
    except _RECOVERABLE_RUR_ERRORS:
        return None
    if result is None or str(result.status).lower() == "unknown":
        return None
    if result.partial and not result.assignments:
        return None
    return result


__all__ = ["finite_points_formula", "try_rur_formula"]
