"""Certified total variation for real univariate functions."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..results import TruthValue
from .common import normalize_domain_set


@dataclass(frozen=True)
class VariationResult:
    bounded_variation: TruthValue
    total_variation: sp.Expr | None
    domain: object
    source: str | None = None
    evidence: tuple[object, ...] = ()


def variation(expr, var, *, domain=S.Reals, return_result=False):
    """Return certified total variation; ``oo`` denotes infinite variation."""
    expr, var = sp.sympify(expr), sp.sympify(var)
    try:
        d = normalize_domain_set(domain, var)
    except (TypeError, ValueError, NotImplementedError):
        d = None
    if not expr.has(var):
        result = VariationResult(TruthValue.TRUE, S.Zero, domain, "constant")
    elif (
        isinstance(d, sp.Interval)
        and d.start.is_finite is True
        and d.end.is_finite is True
    ):
        derivative = sp.diff(expr, var)
        try:
            nonnegative = sp.calculus.util.function_range(
                derivative, var, d
            ).inf.is_nonnegative
            nonpositive = sp.calculus.util.function_range(
                derivative, var, d
            ).sup.is_nonpositive
        except (TypeError, ValueError, NotImplementedError):
            nonnegative = nonpositive = None
        if nonnegative is True or nonpositive is True:
            left = (
                sp.limit(expr, var, d.start, dir="+")
                if d.left_open
                else expr.subs(var, d.start)
            )
            right = (
                sp.limit(expr, var, d.end, dir="-")
                if d.right_open
                else expr.subs(var, d.end)
            )
            result = VariationResult(
                TruthValue.TRUE,
                sp.simplify(sp.Abs(right - left)),
                domain,
                "monotone-interval",
            )
        else:
            try:
                value = sp.integrate(sp.Abs(derivative), (var, d.start, d.end))
            except (TypeError, ValueError, NotImplementedError):
                value = sp.Integral(sp.Abs(derivative), (var, d.start, d.end))
            result = (
                VariationResult(
                    TruthValue.TRUE,
                    sp.simplify(value),
                    domain,
                    "absolute-derivative-integral",
                )
                if not isinstance(value, sp.Integral) and value.is_finite is True
                else VariationResult(TruthValue.UNKNOWN, None, domain, "unresolved")
            )
    else:
        result = VariationResult(TruthValue.UNKNOWN, None, domain, "unresolved")
    return result if return_result else result.total_variation
