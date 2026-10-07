from __future__ import annotations

import contextlib

import mpmath as mp
import sympy as sp

from ..exactness import certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS


def _exact_sign(expr):
    """Return -1, 0, or 1 when SymPy establishes the real sign."""
    expr = sp.sympify(expr)
    with contextlib.suppress(*SYMBOLIC_QUERY_ERRORS):
        expr = sp.simplify(expr)
    if certify_zero(expr) is True:
        return 0
    if expr.is_positive is True:
        return 1
    if expr.is_negative is True:
        return -1
    with contextlib.suppress(*SYMBOLIC_QUERY_ERRORS):
        q = sp.sign(expr)
        if q in (-1, 0, 1):
            return int(q)
    return None


def _iv_rational(q):
    if isinstance(q, sp.Integer):
        return mp.iv.mpf([int(q), int(q)])
    if isinstance(q, sp.Rational):
        return mp.iv.mpf([int(q.p), int(q.p)]) / mp.iv.mpf([int(q.q), int(q.q)])
    raise ValueError


def _interval_eval(expr):
    """Evaluate a numeric exp/log elementary expression with interval arithmetic.

    This intentionally supports a small audited vocabulary. Unsupported heads raise
    ``ValueError`` so callers remain conservative
    to ordinary floating-point evaluation.
    """
    expr = sp.sympify(expr)
    if isinstance(expr, (sp.Integer, sp.Rational)):
        return _iv_rational(expr)
    if expr == sp.pi:
        return mp.iv.pi
    if expr == sp.E:
        return mp.iv.e
    if expr.is_Add:
        out = mp.iv.mpf([0, 0])
        for arg in expr.args:
            out += _interval_eval(arg)
        return out
    if expr.is_Mul:
        out = mp.iv.mpf([1, 1])
        for arg in expr.args:
            out *= _interval_eval(arg)
        return out
    if expr.is_Pow:
        base, exponent = expr.args
        b = _interval_eval(base)
        if exponent.is_Integer:
            return b ** int(exponent)
        if exponent.is_Rational:
            # Only use real positive-base rational powers in the robust sign path.
            if float(b.a) <= 0:
                raise ValueError
            return mp.iv.exp(_iv_rational(exponent) * mp.iv.log(b))
        raise ValueError
    if expr.func is sp.exp:
        return mp.iv.exp(_interval_eval(expr.args[0]))
    if expr.func is sp.log:
        x = _interval_eval(expr.args[0])
        if float(x.a) <= 0:
            raise ValueError
        return mp.iv.log(x)
    if expr.func is sp.sinh:
        x = _interval_eval(expr.args[0])
        return (mp.iv.exp(x) - mp.iv.exp(-x)) / 2
    if expr.func is sp.cosh:
        x = _interval_eval(expr.args[0])
        return (mp.iv.exp(x) + mp.iv.exp(-x)) / 2
    if expr.func is sp.tanh:
        x = _interval_eval(expr.args[0])
        ex = mp.iv.exp(2 * x)
        return (ex - 1) / (ex + 1)
    raise ValueError


def robust_real_sign(expr, *, interval_fallback=True):
    """Return a certified sign for a real numeric expression when possible.

    The exact symbolic route is tried first.  For unresolved numeric expressions in
    an audited exp/log vocabulary, interval arithmetic is used so very large or very
    small exponentials can be compared without overflow/underflow.  A numerical
    interval containing zero is *not* converted into a sign.
    """
    expr = sp.sympify(expr)
    sign = _exact_sign(expr)
    if sign is not None or not interval_fallback:
        return sign
    if expr.free_symbols or expr.is_real is False:
        return None
    try:
        value = _interval_eval(expr)
        if float(value.a) > 0:
            return 1
        if float(value.b) < 0:
            return -1
    except (ValueError, TypeError, OverflowError, ZeroDivisionError):
        return None
    return None


__all__ = ["robust_real_sign"]
