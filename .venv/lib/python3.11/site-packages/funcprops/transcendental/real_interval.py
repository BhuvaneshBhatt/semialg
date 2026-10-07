from __future__ import annotations

import math
from dataclasses import dataclass

import mpmath as mp
import sympy as sp


@dataclass(frozen=True)
class RealIntervalEvaluation:
    expression: sp.Basic
    variable: sp.Symbol
    domain: sp.Interval
    lower: str
    upper: str
    excludes_zero: bool
    strict_sign: int
    source: str = "mpmath-outward-rounded-real-interval"


def _const_iv(expr):
    if isinstance(expr, sp.Integer):
        n = int(expr)
        return mp.iv.mpf([n, n])
    if isinstance(expr, sp.Rational):
        return mp.iv.mpf([int(expr.p), int(expr.p)]) / mp.iv.mpf(
            [int(expr.q), int(expr.q)]
        )
    if expr == sp.pi:
        return mp.iv.pi
    if expr == sp.E:
        return mp.iv.e
    raise ValueError


def _eval(expr, var, x):
    expr = sp.sympify(expr)
    if expr == var:
        return x
    if not expr.has(var) and (
        isinstance(expr, (sp.Integer, sp.Rational)) or expr in (sp.pi, sp.E)
    ):
        return _const_iv(expr)
    if expr.is_Add:
        out = mp.iv.mpf([0, 0])
        for arg in expr.args:
            out += _eval(arg, var, x)
        return out
    if expr.is_Mul:
        out = mp.iv.mpf([1, 1])
        for arg in expr.args:
            out *= _eval(arg, var, x)
        return out
    if expr.is_Pow:
        base, exponent = expr.args
        b = _eval(base, var, x)
        if exponent.is_Integer:
            return b ** int(exponent)
        if exponent.is_Rational and float(b.a) > 0:
            q = _const_iv(exponent)
            return mp.iv.exp(q * mp.iv.log(b))
        raise ValueError
    if expr.func is sp.exp:
        return mp.iv.exp(_eval(expr.args[0], var, x))
    if expr.func is sp.log:
        a = _eval(expr.args[0], var, x)
        if float(a.a) <= 0:
            raise ValueError
        return mp.iv.log(a)
    if expr.func is sp.sin:
        return mp.iv.sin(_eval(expr.args[0], var, x))
    if expr.func is sp.cos:
        return mp.iv.cos(_eval(expr.args[0], var, x))
    if expr.func is sp.tan:
        return mp.iv.tan(_eval(expr.args[0], var, x))
    if expr.func is sp.sinh:
        a = _eval(expr.args[0], var, x)
        return (mp.iv.exp(a) - mp.iv.exp(-a)) / 2
    if expr.func is sp.cosh:
        a = _eval(expr.args[0], var, x)
        return (mp.iv.exp(a) + mp.iv.exp(-a)) / 2
    if expr.func is sp.tanh:
        a = _eval(expr.args[0], var, x)
        e = mp.iv.exp(2 * a)
        return (e - 1) / (e + 1)
    raise ValueError


def certified_real_interval_evaluate(expr, var, domain):
    """Enclose ``expr(domain)`` for a finite exact real interval.

    The routine intentionally supports a bounded elementary vocabulary and returns
    ``None`` outside it.  Its purpose is a conservative root-isolation fallback:
    excluding zero proves a cell root-free, and excluding zero from a derivative
    enclosure proves strict monotonicity on that cell.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if (
        not isinstance(domain, sp.Interval)
        or domain.start in (-sp.oo, sp.oo)
        or domain.end in (-sp.oo, sp.oo)
        or domain.start.is_real is False
        or domain.end.is_real is False
    ):
        return None
    try:
        a = _const_iv(domain.start)
        b = _const_iv(domain.end)
        x = mp.iv.mpf([a.a, b.b])
        value = _eval(expr, var, x)
        lo = math.nextafter(float(mp.mpf(value.a._mpi_[0])), -math.inf)
        hi = math.nextafter(float(mp.mpf(value.b._mpi_[1])), math.inf)
    except (ValueError, TypeError, OverflowError, ZeroDivisionError):
        return None
    sign = 1 if lo > 0 else -1 if hi < 0 else 0
    return RealIntervalEvaluation(
        expr,
        var,
        domain,
        repr(lo),
        repr(hi),
        sign != 0,
        sign,
    )


__all__ = ["RealIntervalEvaluation", "certified_real_interval_evaluate"]
