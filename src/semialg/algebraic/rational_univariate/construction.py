from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp
from algroots.errors import RationalUnivariateError as _AlgrootsRationalUnivariateError
from algroots.rational_univariate import compute_rational_univariate_representation

from ..cache import CACHE, expr_key
from ..groebner_utils import _backend_policy
from .representation import RationalUnivariateError, RationalUnivariateRepresentation


@_backend_policy
def compute_rur(
    polynomials: Iterable[sp.Expr],
    variables: Sequence[sp.Symbol],
    parameter: sp.Symbol | None = None,
    *,
    max_separating_attempts: int = 64,
    domain=None,
    groebner_backend="native",
    groebner_limits=None,
) -> RationalUnivariateRepresentation:
    """Compute a RUR using algroots as the canonical equality-only backend.

    ``domain`` declares an exact characteristic-zero coefficient field, including
    rational parameter fields. ``groebner_backend`` and ``groebner_limits`` select verified construction
    and participate in the cache key. Semialg retains only its cache and semialgebraic integration around the
    canonical algroots finite-quotient implementation.
    """

    from algroots.basis import _validate_backend

    _validate_backend(groebner_backend, groebner_limits)
    if (
        isinstance(max_separating_attempts, bool)
        or not isinstance(max_separating_attempts, int)
        or max_separating_attempts < 1
    ):
        raise ValueError("max_separating_attempts must be a positive integer")
    from algroots.domains import coefficient_domain

    domain = coefficient_domain(domain) if domain is not None else None
    variable_tuple = tuple(variables)
    raw_polynomials = tuple(sp.sympify(poly) for poly in polynomials)
    cache_key = (
        tuple(expr_key(poly) for poly in raw_polynomials),
        variable_tuple,
        parameter,
        int(max_separating_attempts),
        "algroots",
        domain,
        groebner_backend,
        groebner_limits,
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
            domain=domain,
            groebner_backend=groebner_backend,
            groebner_limits=groebner_limits,
        )
    except _AlgrootsRationalUnivariateError as exc:
        raise RationalUnivariateError(str(exc)) from exc
    CACHE.rur.put(cache_key, result)
    return result


__all__ = ["compute_rur"]
