"""Structured integrability classification without conflating integral notions."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..results import Continuity, TruthValue
from ..singularities import continuous
from .common import normalize_domain_set


@dataclass(frozen=True)
class IntegrabilityResult:
    riemann: TruthValue = TruthValue.UNKNOWN
    locally_integrable: TruthValue = TruthValue.UNKNOWN
    absolutely_integrable: TruthValue = TruthValue.UNKNOWN
    improper_integrable: TruthValue = TruthValue.UNKNOWN
    integral: sp.Expr | None = None
    absolute_integral: sp.Expr | None = None
    domain: object = S.Reals
    exceptional_loci: tuple[object, ...] = ()
    evidence: tuple[object, ...] = ()


def _convergence(integral):
    value = _evaluated(integral)
    if value is None:
        return TruthValue.UNKNOWN
    if value.has(S.Infinity, S.NegativeInfinity, sp.zoo, sp.nan):
        return TruthValue.FALSE
    return TruthValue.TRUE if value.is_finite is True else TruthValue.UNKNOWN


def _evaluated(integral):
    try:
        value = integral.doit()
    except (TypeError, ValueError, NotImplementedError):
        return None
    return None if isinstance(value, sp.Integral) else value


def integrability(expr, var, *, domain=S.Reals, return_result=False):
    """Classify Riemann, local, absolute (L1), and improper integrability."""
    expr, var = sp.sympify(expr), sp.sympify(var)
    try:
        d = normalize_domain_set(domain, var)
    except (TypeError, ValueError, NotImplementedError):
        d = None
    riemann = local = absolute = improper = TruthValue.UNKNOWN
    integral_value = absolute_value = None
    evidence = []
    if isinstance(d, sp.Interval):
        finite = d.start.is_finite is True and d.end.is_finite is True
        continuity = continuous(expr, var, domain=d)
        if (
            finite
            and not d.left_open
            and not d.right_open
            and continuity is Continuity.CONTINUOUS
        ):
            riemann = local = absolute = TruthValue.TRUE
            evidence.append("continuous-on-compact-interval")
        elif continuity is Continuity.CONTINUOUS:
            local = TruthValue.TRUE
            evidence.append("continuous-on-domain")
        integral = sp.Integral(expr, (var, d.start, d.end))
        absolute_integral = sp.Integral(sp.Abs(expr), (var, d.start, d.end))
        if not finite or d.left_open or d.right_open:
            improper = _convergence(integral)
        absolute = _convergence(absolute_integral)
        integral_value = (
            _evaluated(integral)
            if improper is TruthValue.TRUE or riemann is TruthValue.TRUE
            else None
        )
        absolute_value = (
            _evaluated(absolute_integral) if absolute is TruthValue.TRUE else None
        )
    elif d is S.Reals:
        if continuous(expr, var, domain=S.Reals) is Continuity.CONTINUOUS:
            local = TruthValue.TRUE
        integral = sp.Integral(expr, (var, -S.Infinity, S.Infinity))
        absolute_integral = sp.Integral(sp.Abs(expr), (var, -S.Infinity, S.Infinity))
        improper, absolute = _convergence(integral), _convergence(absolute_integral)
        integral_value = _evaluated(integral) if improper is TruthValue.TRUE else None
        absolute_value = (
            _evaluated(absolute_integral) if absolute is TruthValue.TRUE else None
        )
    result = IntegrabilityResult(
        riemann,
        local,
        absolute,
        improper,
        integral_value,
        absolute_value,
        domain,
        evidence=tuple(evidence),
    )
    return result
