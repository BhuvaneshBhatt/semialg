from __future__ import annotations

import sympy as sp

from ._common import RECOVERABLE_ERRORS


def linear_equality_candidate(
    expanded_difference: sp.Expr,
    variable: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr] | None:
    """Return ``(coefficient, numerator, denominator, replacement)`` for a
    linear equality in ``variable``.

    ``expanded_difference`` is assumed to be the already-expanded ``lhs-rhs``;
    accepting that value avoids repeating expansion in callers that scan several
    variables of the same equality.
    """

    try:
        poly = sp.Poly(expanded_difference, variable)
    except RECOVERABLE_ERRORS:
        return None
    if poly.degree() != 1:
        return None
    coefficient = sp.expand(poly.coeff_monomial(variable))
    constant = sp.expand(poly.coeff_monomial(1))
    if coefficient == 0 or coefficient.has(variable) or constant.has(variable):
        return None
    numerator = sp.expand(-constant)
    denominator = coefficient
    replacement = sp.simplify(numerator / denominator)
    return coefficient, numerator, denominator, replacement


__all__ = ["linear_equality_candidate"]
