from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp
from algroots.errors import RationalUnivariateError as _AlgrootsRationalUnivariateError
from algroots.rational_univariate import solve_rur_points as _algroots_solve_rur_points
from algroots.rational_univariate import (
    solve_rur_representation as _algroots_solve_rur_representation,
)

from .construction import compute_rur
from .representation import (
    RationalUnivariateError,
    RationalUnivariatePoint,
    RationalUnivariateRepresentation,
)


def solve_rur_representation(
    representation: RationalUnivariateRepresentation,
    *,
    real: bool = True,
) -> tuple[tuple[sp.Expr, ...], ...]:
    """Solve an existing RUR through the canonical algroots backend."""

    try:
        return _algroots_solve_rur_representation(representation, real=real)
    except _AlgrootsRationalUnivariateError as exc:
        raise RationalUnivariateError(str(exc)) from exc


def solve_with_rur(
    polynomials: Iterable[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    real: bool = True,
    parameter: sp.Symbol | None = None,
    max_separating_attempts: int = 64,
) -> tuple[tuple[sp.Expr, ...], ...]:
    """Return distinct exact solutions obtained from an algroots RUR."""

    representation = compute_rur(
        polynomials,
        variables,
        parameter,
        max_separating_attempts=max_separating_attempts,
    )
    return solve_rur_representation(representation, real=real)


def solve_rur_points(
    representation: RationalUnivariateRepresentation,
    *,
    real: bool = True,
) -> tuple[RationalUnivariatePoint, ...]:
    """Return Semialg-enriched RUR points backed by algroots root isolation."""

    try:
        points = _algroots_solve_rur_points(representation, real=real)
    except _AlgrootsRationalUnivariateError as exc:
        raise RationalUnivariateError(str(exc)) from exc
    return tuple(RationalUnivariatePoint(point.representation, point.root) for point in points)


__all__ = ["solve_with_rur", "solve_rur_representation", "solve_rur_points"]
