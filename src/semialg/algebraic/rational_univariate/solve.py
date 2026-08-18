from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp

from .construction import compute_rational_univariate_representation
from .representation import (
    RationalUnivariatePoint,
    RationalUnivariateRepresentation,
)


def solve_rur_representation(
    representation: RationalUnivariateRepresentation,
    *,
    real: bool = True,
) -> tuple[tuple[sp.Expr, ...], ...]:
    """Return distinct exact solutions from an existing RUR."""

    if representation.defining_polynomial.degree() <= 0:
        return tuple()

    t = representation.parameter
    defining = representation.defining_polynomial
    coordinate_polys = representation.normalized_coordinate_polynomials()
    if real:
        roots = tuple(sp.real_roots(defining.as_expr()))
    else:
        roots = tuple(sp.Poly(defining.as_expr(), t, domain=sp.QQ).all_roots())
    solutions: list[tuple[sp.Expr, ...]] = []
    seen_keys: set[tuple[str, ...]] = set()
    for root in roots:
        point = tuple(sp.cancel(poly.as_expr().subs(t, root)) for poly in coordinate_polys)
        key = tuple(sp.sstr(coord) for coord in point)
        if key in seen_keys:
            continue
        if any(
            all(sp.simplify(a - b) == 0 for a, b in zip(point, old, strict=True))
            for old in solutions
        ):
            continue
        seen_keys.add(key)
        solutions.append(point)
    return tuple(solutions)


def solve_zero_dimensional_system_with_rur(
    polynomials: Iterable[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    real: bool = True,
    parameter: sp.Symbol | None = None,
    max_separating_attempts: int = 64,
) -> tuple[tuple[sp.Expr, ...], ...]:
    """Return distinct exact solutions obtained from a RUR representation."""

    representation = compute_rational_univariate_representation(
        polynomials, variables, parameter, max_separating_attempts=max_separating_attempts
    )
    return solve_rur_representation(representation, real=real)


def solve_rur_points(
    representation: RationalUnivariateRepresentation,
    *,
    real: bool = True,
) -> tuple[RationalUnivariatePoint, ...]:
    """Return distinct solutions as RUR parameter-root points."""

    if representation.defining_polynomial.degree() <= 0:
        return tuple()
    t = representation.parameter
    defining = representation.defining_polynomial
    roots = (
        tuple(sp.real_roots(defining.as_expr()))
        if real
        else tuple(sp.Poly(defining.as_expr(), t, domain=sp.QQ).all_roots())
    )
    points: list[RationalUnivariatePoint] = []
    seen: set[tuple[str, ...]] = set()
    for root in roots:
        point = RationalUnivariatePoint(representation, root)
        key = tuple(sp.sstr(coord) for coord in point.coordinates)
        if key in seen:
            continue
        seen.add(key)
        points.append(point)
    return tuple(points)
