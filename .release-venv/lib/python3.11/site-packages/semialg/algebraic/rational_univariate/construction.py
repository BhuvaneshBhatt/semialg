from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp
from algroots.errors import RationalUnivariateError as _AlgrootsRationalUnivariateError
from algroots.rational_univariate import compute_rational_univariate_representation

from ..cache import CACHE, expr_key
from .representation import RationalUnivariateError, RationalUnivariateRepresentation


def compute_rur(
    polynomials: Iterable[sp.Expr],
    variables: Sequence[sp.Symbol],
    parameter: sp.Symbol | None = None,
    *,
    max_separating_attempts: int = 64,
) -> RationalUnivariateRepresentation:
    """Compute a RUR using algroots as the canonical equality-only backend.

    Semialg retains only its cache and semialgebraic integration around the
    canonical algroots finite-quotient implementation.
    """

    variable_tuple = tuple(variables)
    raw_polynomials = tuple(sp.sympify(poly) for poly in polynomials)
    cache_key = (
        tuple(expr_key(poly) for poly in raw_polynomials),
        variable_tuple,
        parameter,
        int(max_separating_attempts),
        "algroots",
    )
    cached = CACHE.rur.get(cache_key)
    if cached is not None:
        CACHE.stats.rur_hits += 1
        return cached  # type: ignore[return-value]
    CACHE.stats.rur_misses += 1
    try:
        result = compute_rational_univariate_representation(
            raw_polynomials,
            variable_tuple,
            parameter,
            max_separating_attempts=max_separating_attempts,
        )
    except _AlgrootsRationalUnivariateError as exc:
        raise RationalUnivariateError(str(exc)) from exc
    CACHE.rur.put(cache_key, result)
    return result


__all__ = ["compute_rur"]
