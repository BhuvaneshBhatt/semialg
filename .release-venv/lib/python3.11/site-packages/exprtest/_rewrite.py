"""Small size-bounded exact rewrite registry for the fast oracle."""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from . import _config as cfg
from ._fast import quick_reduce


def _one(term: sp.Expr) -> sp.Expr:
    if term.func is sp.Abs and len(term.args) == 1:
        arg = term.args[0]
        if arg.is_nonnegative is True:
            return arg
        if arg.is_nonpositive is True:
            return -arg
    if term.func is sp.sign and len(term.args) == 1:
        arg = term.args[0]
        if arg.is_positive is True:
            return sp.Integer(1)
        if arg.is_negative is True:
            return sp.Integer(-1)
        if arg.is_zero is True:
            return sp.Integer(0)
    if term.func is sp.conjugate and len(term.args) == 1:
        arg = term.args[0]
        if arg.is_real is True:
            return arg
        if arg.is_imaginary is True:
            return -arg
    if term.func is sp.re and len(term.args) == 1:
        arg = term.args[0]
        if arg.is_real is True:
            return arg
        if arg.is_imaginary is True:
            return sp.Integer(0)
    if term.func is sp.im and len(term.args) == 1:
        arg = term.args[0]
        if arg.is_real is True:
            return sp.Integer(0)
        if arg.is_imaginary is True:
            return -sp.I * arg
    if term.func is sp.log and len(term.args) == 1 and term.args[0] == 1:
        return sp.Integer(0)
    if term.func is sp.exp and len(term.args) == 1 and term.args[0] == 0:
        return sp.Integer(1)
    if term.func in (sp.sinh, sp.tanh) and len(term.args) == 1 and term.args[0] == 0:
        return sp.Integer(0)
    if term.func is sp.cosh and len(term.args) == 1 and term.args[0] == 0:
        return sp.Integer(1)
    if term.is_Pow:
        base, exponent = term.args
        if base in (sp.Integer(1), sp.Integer(-1), sp.I, -sp.I) and exponent.is_Integer:
            return base ** (
                int(exponent) % 4 if base in (sp.I, -sp.I) else int(exponent) % 2
            )
    return term


@lru_cache(maxsize=cfg.EXACT_REWRITE_CACHE_SIZE)
def exact_rewrite(term: sp.Expr) -> sp.Expr:
    """Apply exact local rewrites while keeping growth tightly bounded."""
    term = sp.sympify(term)
    if term.count_ops() > cfg.EXACT_REWRITE_MAX_OPS:
        return term
    if not term.args:
        return term
    args = tuple(exact_rewrite(arg) for arg in term.args)
    try:
        rebuilt = term.func(*args, evaluate=False)
    except (TypeError, ValueError):
        try:
            rebuilt = term.func(*args)
        except (TypeError, ValueError):
            rebuilt = term
    reduced = _one(rebuilt)
    if reduced.count_ops() > rebuilt.count_ops() + 1:
        return rebuilt
    return quick_reduce(reduced)


def clear_rewrite_cache() -> None:
    exact_rewrite.cache_clear()
