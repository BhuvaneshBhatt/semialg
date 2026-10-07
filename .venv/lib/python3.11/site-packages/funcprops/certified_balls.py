"""Certified real-ball facade over outward interval evaluation."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .transcendental.real_interval import certified_real_interval_evaluate


@dataclass(frozen=True)
class RealBallCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    domain: sp.Interval
    lower: sp.Expr
    upper: sp.Expr
    backend: str = "outward-interval"
    precision: int | None = None

    @property
    def excludes_zero(self):
        return self.lower.is_positive is True or self.upper.is_negative is True

    @property
    def sign(self):
        return (
            1
            if self.lower.is_positive is True
            else -1
            if self.upper.is_negative is True
            else 0
        )

    def verify(self):
        r = certified_real_interval_evaluate(
            self.expression, self.variable, self.domain
        )
        if r is None:
            return False
        try:
            lo = sp.Float(r.lower)
            hi = sp.Float(r.upper)
        except (TypeError, ValueError):
            return False
        return lo >= self.lower and hi <= self.upper


def certified_real_ball(expr, var, domain, *, precision=None):
    r = certified_real_interval_evaluate(expr, var, domain)
    if r is None:
        return None
    try:
        lo = sp.Float(r.lower)
        hi = sp.Float(r.upper)
    except (TypeError, ValueError):
        return None
    return RealBallCertificate(
        sp.sympify(expr), sp.sympify(var), domain, lo, hi, precision=precision
    )
