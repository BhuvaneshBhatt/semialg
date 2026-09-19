"""Cheap exact special-function values for the fast oracle."""

from __future__ import annotations

from functools import lru_cache
from math import comb as py_comb
from math import factorial as py_factorial

import sympy as sp

from . import _config as cfg
from ._cost import stage_allowed
from ._errors import EXACT_METHOD_ERRORS


def _small_int(arg: sp.Expr, *, nonnegative: bool = False):
    if not arg.is_Integer:
        return None
    value = int(arg)
    if nonnegative and value < 0:
        return None
    if abs(value) > cfg.EXACT_SPECIAL_MAX_INDEX:
        return None
    return value


def _gamma_value(arg: sp.Expr):
    n = _small_int(arg)
    if n is not None and n > 0:
        return sp.Integer(py_factorial(n - 1))
    if arg.is_Rational and arg.q == 2:
        # Gamma(1/2 + n), including negative half-integers away from poles.
        shift = arg - sp.Rational(1, 2)
        n = _small_int(shift)
        if n is None:
            return None
        if n >= 0:
            return sp.Rational(py_factorial(2 * n), 4**n * py_factorial(n)) * sp.sqrt(
                sp.pi
            )
        k = -n
        return sp.Rational((-4) ** k * py_factorial(k), py_factorial(2 * k)) * sp.sqrt(
            sp.pi
        )
    return None


def _zeta_value(arg: sp.Expr):
    n = _small_int(arg)
    if n is None or n == 1:
        return None
    if n <= 0:
        return -sp.bernoulli(1 - n) / sp.Integer(1 - n)
    if n % 2 == 0:
        k = n // 2
        return (
            (-1) ** (k + 1)
            * sp.bernoulli(2 * k)
            * (2 * sp.pi) ** (2 * k)
            / (2 * sp.Integer(py_factorial(2 * k)))
        )
    return None


def _factorial2_value(arg: sp.Expr):
    n = _small_int(arg)
    if n is None or n < -1:
        return None
    if n in (-1, 0, 1):
        return sp.Integer(1)
    value = 1
    for k in range(n, 0, -2):
        value *= k
    return sp.Integer(value)


def _binomial_value(n: sp.Expr, k: sp.Expr):
    ni = _small_int(n)
    ki = _small_int(k, nonnegative=True)
    if ni is None or ki is None:
        return None
    if ni >= 0:
        return sp.Integer(0 if ki > ni else py_comb(ni, ki))
    # Generalized integer binomial: C(n,k)=(-1)^k C(k-n-1,k).
    return sp.Integer((-1) ** ki * py_comb(ki - ni - 1, ki))


def _rising_value(base: sp.Expr, count: sp.Expr):
    n = _small_int(count, nonnegative=True)
    if n is None or not base.is_Rational:
        return None
    value = sp.Integer(1)
    for k in range(n):
        value *= base + k
    return value


def _falling_value(base: sp.Expr, count: sp.Expr):
    n = _small_int(count, nonnegative=True)
    if n is None or not base.is_Rational:
        return None
    value = sp.Integer(1)
    for k in range(n):
        value *= base - k
    return value


def _harmonic_value(args: tuple[sp.Expr, ...]):
    if not args or len(args) > 2:
        return None
    n = _small_int(args[0], nonnegative=True)
    if n is None:
        return None
    order = sp.Integer(1) if len(args) == 1 else args[1]
    m = _small_int(order, nonnegative=True)
    if m is None or m == 0:
        return sp.Integer(n)
    return sum((sp.Rational(1, k**m) for k in range(1, n + 1)), sp.Integer(0))


def _one_value(term: sp.Expr):
    if term.func is sp.gamma and len(term.args) == 1:
        return _gamma_value(term.args[0])
    if term.func is sp.zeta and len(term.args) == 1:
        return _zeta_value(term.args[0])
    if term.func is sp.factorial and len(term.args) == 1:
        n = _small_int(term.args[0], nonnegative=True)
        return None if n is None else sp.Integer(py_factorial(n))
    if term.func is sp.factorial2 and len(term.args) == 1:
        return _factorial2_value(term.args[0])
    if term.func is sp.binomial and len(term.args) == 2:
        return _binomial_value(*term.args)
    if term.func is sp.rf and len(term.args) == 2:
        return _rising_value(*term.args)
    if term.func is sp.ff and len(term.args) == 2:
        return _falling_value(*term.args)
    if term.func is sp.harmonic:
        return _harmonic_value(term.args)
    if term.func is sp.beta and len(term.args) == 2:
        a, b = term.args
        ga = _gamma_value(a)
        gb = _gamma_value(b)
        gab = _gamma_value(a + b)
        if ga is not None and gb is not None and gab is not None and gab != 0:
            return ga * gb / gab
    return None


@lru_cache(maxsize=2048)
def _special_cached(term: sp.Expr) -> sp.Expr:
    if not term.args:
        return term
    args = tuple(_special_cached(arg) for arg in term.args)
    try:
        rebuilt = term.func(*args, evaluate=False)
    except EXACT_METHOD_ERRORS:
        try:
            rebuilt = term.func(*args)
        except EXACT_METHOD_ERRORS:
            rebuilt = term
    value = _one_value(rebuilt)
    return rebuilt if value is None else value


def reduce_special_values(term: sp.Expr) -> sp.Expr:
    """Rewrite a bounded registry of rigorously known exact special values."""
    term = sp.sympify(term)
    if not stage_allowed(term, "special"):
        return term
    return _special_cached(term)
