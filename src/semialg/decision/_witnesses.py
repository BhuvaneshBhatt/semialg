from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from ..instances.real_fallbacks import satisfies_formula
from ..sampling import sample_points

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    RuntimeError,
    sp.PolynomialError,
)


def validate_witness(
    formula: sp.Expr,
    point: Mapping[sp.Symbol, sp.Expr] | None,
    variables: Sequence[sp.Symbol],
) -> Mapping[sp.Symbol, sp.Expr] | None:
    if point is None:
        return None
    variables = tuple(variables)
    normalized = {var: sp.sympify(point[var]) for var in variables if var in point}
    if len(normalized) != len(variables):
        return None
    try:
        if satisfies_formula(formula, normalized):
            return normalized
    except _RECOVERABLE_ERRORS:
        pass
    return None


def find_validated_witness(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    strategy: str | None = None,
) -> Mapping[sp.Symbol, sp.Expr] | None:
    try:
        points = sample_points(formula, variables, count=1, strategy=strategy or "auto")
    except _RECOVERABLE_ERRORS:
        points = ()
    for point in points:
        validated = validate_witness(formula, point, variables)
        if validated is not None:
            return validated
    return None


__all__ = ["find_validated_witness", "validate_witness"]
