"""Exact Hurwitz stability regions for real polynomial characteristic equations."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..normalization import normalize_variables
from ..symbol_resolution import resolve_symbol


@dataclass(frozen=True)
class PolynomialStabilityResult:
    """Strict left-half-plane stability conditions for a real polynomial."""

    polynomial: sp.Expr
    variable: sp.Symbol
    parameters: tuple[sp.Symbol, ...]
    degree: int
    condition: sp.Expr
    hurwitz_matrix: sp.Matrix
    determinants: tuple[sp.Expr, ...]
    method: str = "routh_hurwitz"
    certified: bool = True
    diagnostics: Mapping[str, object] = field(default_factory=dict)


def _hurwitz_matrix(coeffs: Sequence[sp.Expr]) -> sp.Matrix:
    """Build the square Hurwitz matrix from descending polynomial coefficients."""

    degree = len(coeffs) - 1

    def coefficient(power: int) -> sp.Expr:
        if power < 0 or power > degree:
            return sp.S.Zero
        return coeffs[degree - power]

    rows: list[list[sp.Expr]] = []
    for row in range(degree):
        offset = row // 2
        odd_row = row % 2 == 0
        values: list[sp.Expr] = []
        for col in range(degree):
            if odd_row:
                power = degree - 1 - 2 * (col - offset)
            else:
                power = degree - 2 * (col - offset)
            values.append(coefficient(power))
        rows.append(values)
    return sp.Matrix(rows)


def _hurwitz_data(coeffs: Sequence[sp.Expr]) -> tuple[sp.Matrix, tuple[sp.Expr, ...]]:
    matrix = _hurwitz_matrix(coeffs)
    dets = tuple(sp.factor(matrix[:size, :size].det()) for size in range(1, matrix.rows + 1))
    return matrix, dets


def polynomial_stability_analysis(
    polynomial: sp.Expr,
    variable: sp.Symbol | str,
    parameters: Sequence[sp.Symbol | str] | None = None,
) -> PolynomialStabilityResult:
    """Return exact strict Hurwitz-stability conditions for a real polynomial.

    The returned condition characterizes parameter values for which every root
    of ``polynomial`` has strictly negative real part.  Parameter values for
    which the leading coefficient vanishes are excluded because the polynomial
    degree changes there.
    """

    expr = sp.expand(sp.sympify(polynomial))
    var = resolve_symbol(variable, context=(expr,))
    try:
        poly = sp.Poly(expr, var)
    except sp.PolynomialError as exc:
        raise ValueError(
            "characteristic equation must be polynomial in the stability variable"
        ) from exc
    degree = int(poly.degree())
    if degree < 1:
        raise ValueError("characteristic polynomial must have positive degree")
    coeffs = tuple(sp.sympify(item) for item in poly.all_coeffs())
    if any(item.is_real is False for item in coeffs):
        raise ValueError("Hurwitz stability requires real polynomial coefficients")

    free = tuple(sorted(expr.free_symbols - {var}, key=lambda sym: (sym.name, sp.srepr(sym))))
    if parameters is None:
        params = free
    else:
        params = normalize_variables(parameters, expr, append_context_symbols=False)
        missing = set(free) - set(params)
        if missing:
            names = ", ".join(sorted(sym.name for sym in missing))
            raise ValueError(f"parameters do not include all coefficient symbols: {names}")

    matrix, dets = _hurwitz_data(coeffs)
    lead = coeffs[0]
    pos_condition = sp.And(lead > 0, *(det > 0 for det in dets))

    neg_coeffs = tuple(-item for item in coeffs)
    _, neg_dets = _hurwitz_data(neg_coeffs)
    neg_condition = sp.And(lead < 0, *(det > 0 for det in neg_dets))

    if lead.is_positive is True:
        condition = sp.And(*(det > 0 for det in dets))
    elif lead.is_negative is True:
        condition = sp.And(*(det > 0 for det in neg_dets))
    else:
        condition = sp.Or(pos_condition, neg_condition)
    condition = sp.simplify_logic(condition)

    return PolynomialStabilityResult(
        polynomial=expr,
        variable=var,
        parameters=params,
        degree=degree,
        condition=condition,
        hurwitz_matrix=matrix,
        determinants=dets,
        diagnostics={"leading_coefficient": lead, "strict_stability": True},
    )


def polynomial_stability_region(
    polynomial: sp.Expr,
    variable: sp.Symbol | str,
    parameters: Sequence[sp.Symbol | str] | None = None,
) -> sp.Expr:
    """Return the exact parameter region for strict Hurwitz stability."""

    return polynomial_stability_analysis(polynomial, variable, parameters).condition


__all__ = [
    "PolynomialStabilityResult",
    "polynomial_stability_analysis",
    "polynomial_stability_region",
]
