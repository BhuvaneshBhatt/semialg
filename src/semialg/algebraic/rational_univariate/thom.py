from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp

from ..thom import ThomEncoding, thom_encoding
from .representation import RationalUnivariatePoint


def rur_parameter_thom_encoding(point: RationalUnivariatePoint) -> ThomEncoding:
    """Return the Thom encoding certifying a RUR point's parameter root."""

    representation = point.representation
    return thom_encoding(representation.defining_polynomial, point.root)


def sign_at_rur_point(poly: sp.Poly | sp.Expr, point: RationalUnivariatePoint) -> int:
    """Determine a polynomial sign at a RUR point through univariate reduction."""

    representation = point.representation
    t = representation.parameter
    coordinates = representation.normalized_coordinate_polynomials()
    coordinate_map = dict(
        zip(representation.variables, (p.as_expr() for p in coordinates), strict=True)
    )
    expr = poly.as_expr() if isinstance(poly, sp.Poly) else sp.sympify(poly)
    numerator, denominator = sp.fraction(sp.cancel(expr.subs(coordinate_map)))
    defining = representation.defining_polynomial
    domain = defining.domain
    try:
        numerator_poly = sp.Poly(numerator, t, domain=domain).rem(defining)
        denominator_poly = sp.Poly(denominator, t, domain=domain).rem(defining)
    except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed, ValueError) as exc:
        raise ValueError(
            "RUR sign determination requires a rational function in the RUR parameter"
        ) from exc
    encoding = rur_parameter_thom_encoding(point)
    numerator_sign = encoding.sign_of(numerator_poly)
    if numerator_sign == 0:
        return 0
    denominator_sign = encoding.sign_of(denominator_poly)
    if denominator_sign == 0:
        raise ValueError("expression is undefined at the represented RUR point")
    return numerator_sign * denominator_sign


def sign_conditions_at_rur_points(
    polynomials: Iterable[sp.Poly | sp.Expr],
    points: Sequence[RationalUnivariatePoint],
) -> tuple[tuple[int, ...], ...]:
    """Return exact sign conditions of ``polynomials`` at each RUR point."""

    polynomial_tuple = tuple(polynomials)
    return tuple(
        tuple(sign_at_rur_point(poly, point) for poly in polynomial_tuple) for point in points
    )


__all__ = ["rur_parameter_thom_encoding", "sign_at_rur_point", "sign_conditions_at_rur_points"]
