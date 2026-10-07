from __future__ import annotations

import sympy as sp

from ..results import (
    ConditionalExpression,
    Continuity,
    Differentiability,
    PropertyResult,
)
from ..stratification import normalize_cases
from .smoothness import smoothness as _smoothness


def _differentiability_implications(
    primary: Differentiability,
) -> tuple[Differentiability, ...]:
    """Return strictly weaker differentiability properties implied by ``primary``."""
    if primary is Differentiability.SMOOTH:
        return (Differentiability.CK, Differentiability.DIFFERENTIABLE)
    if primary is Differentiability.CK:
        return (Differentiability.DIFFERENTIABLE,)
    return ()


def _classify_smoothness(smoothness, order):
    """Map a regularity certificate to the public differentiability class."""
    verified = smoothness.differentiability_order
    if order is not None:
        if verified is sp.oo or verified == sp.oo:
            return Differentiability.SMOOTH
        if isinstance(verified, int) and verified >= order:
            return (
                Differentiability.DIFFERENTIABLE if order == 1 else Differentiability.CK
            )
        if smoothness.continuity is Continuity.DISCONTINUOUS or (
            smoothness.smooth is False
            and isinstance(verified, int)
            and verified < order
        ):
            return Differentiability.NOT_DIFFERENTIABLE
        return Differentiability.UNKNOWN
    if verified is sp.oo or verified == sp.oo or smoothness.smooth is True:
        return Differentiability.SMOOTH
    if smoothness.continuity is Continuity.DISCONTINUOUS:
        return Differentiability.NOT_DIFFERENTIABLE
    if isinstance(verified, int) and verified >= 2:
        return Differentiability.CK
    if verified == 1:
        return Differentiability.DIFFERENTIABLE
    if smoothness.smooth is False and verified == 0:
        return Differentiability.NOT_DIFFERENTIABLE
    return Differentiability.UNKNOWN


def differentiability(
    expr,
    var,
    *,
    domain=True,
    order=None,
    all_properties=False,
    return_result=False,
    max_order=6,
):
    """Classify differentiability on the requested real domain.

    Unlike a syntactic ``diff`` test, this query requires the relevant
    derivatives to be continuous on the whole requested domain. ``order=k``
    asks whether the function is at least :math:`C^k`; without an order the
    strongest certified class is returned. Parameter strata are preserved.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if order is not None:
        if not isinstance(order, int) or isinstance(order, bool) or order < 1:
            raise ValueError("order must be a positive integer or None")
        probe_order = max(max_order, order)
    else:
        probe_order = max_order
    if (
        not isinstance(probe_order, int)
        or isinstance(probe_order, bool)
        or probe_order < 1
    ):
        raise ValueError("max_order must be a positive integer")

    # Differentiability is a derivative-based property.  Populate the shared raw
    # derivative cache before the higher-level smoothness analysis so callers
    # consistently reuse the unfactored derivative proof object.  Keep
    # factorization lazy: several structural smoothness rules (for example Abs)
    # can otherwise bypass the derivative layer entirely.
    from ..proof_cache import cached_derivative

    cached_derivative(expr, var, order=1)
    smoothness = _smoothness(
        expr, var, domain=domain, max_order=probe_order, return_result=True
    )
    if isinstance(smoothness, PropertyResult):
        out = normalize_cases(
            (
                ConditionalExpression(_classify_smoothness(c.value, order), c.condition)
                for c in smoothness.cases
            ),
            source="differentiability-parameter-stratification",
            evidence=smoothness.evidence,
        )
        return out

    primary = _classify_smoothness(smoothness, order)
    if all_properties or return_result:
        return PropertyResult(
            primary=primary,
            all_properties=_differentiability_implications(primary),
            source="smoothness-analysis",
            exceptions=smoothness.exceptional_loci,
            evidence=(smoothness,),
        )
    return primary
