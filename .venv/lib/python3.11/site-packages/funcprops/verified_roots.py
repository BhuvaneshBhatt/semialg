"""Verified real root exclusion/uniqueness by interval signs and derivative enclosure."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .certified_balls import certified_real_ball


@dataclass(frozen=True)
class RealRootInterval:
    expression: sp.Expr
    variable: sp.Symbol
    interval: sp.Interval
    root_count: int | None
    unique: bool
    source: str

    def verify(self):
        f = self.expression
        x = self.variable
        interval = self.interval
        b = certified_real_ball(f, x, interval)
        if self.root_count == 0:
            return b is not None and b.excludes_zero
        if self.root_count == 1:
            db = certified_real_ball(sp.diff(f, x), x, interval)
            if db is None or not db.excludes_zero:
                return False
            fa = sp.simplify(f.subs(x, interval.start))
            fb = sp.simplify(f.subs(x, interval.end))
            return (
                fa.is_zero is True
                or fb.is_zero is True
                or (fa.is_positive is True and fb.is_negative is True)
                or (fa.is_negative is True and fb.is_positive is True)
            )
        return False


def verify_real_root_interval(expr, var, interval):
    f = sp.sympify(expr)
    x = sp.sympify(var)
    b = certified_real_ball(f, x, interval)
    if b is not None and b.excludes_zero:
        return RealRootInterval(f, x, interval, 0, False, "interval-exclusion")
    db = certified_real_ball(sp.diff(f, x), x, interval)
    if db is not None and db.excludes_zero:
        fa = sp.simplify(f.subs(x, interval.start))
        fb = sp.simplify(f.subs(x, interval.end))
        crossing = (
            fa.is_zero is True
            or fb.is_zero is True
            or (fa.is_positive is True and fb.is_negative is True)
            or (fa.is_negative is True and fb.is_positive is True)
        )
        if crossing:
            return RealRootInterval(f, x, interval, 1, True, "monotone-sign-change")
    return RealRootInterval(f, x, interval, None, False, "inconclusive")
