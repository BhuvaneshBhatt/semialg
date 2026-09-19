"""Small theorem-only exponential independence checks."""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from . import _config as cfg
from ._algebraic_bounds import algebraic_gap_test
from ._domains import NumberKind, number_kind
from ._result import Verdict
from ._sqrt_sum import square_root_sum_test
from ._timing import run_with_time_budget


def _alg_coeff(term: sp.Expr) -> bool:
    return number_kind(term) is NumberKind.ALGEBRAIC


@lru_cache(maxsize=cfg.EXACT_BOUND_CACHE_SIZE)
def _algebraic_nonzero(term: sp.Expr) -> bool | None:
    """Prove an algebraic exponent difference nonzero with a cheap-first path."""
    term = sp.sympify(term)
    if term.is_Rational:
        return term != 0
    ops = int(term.count_ops())
    if ops <= 4:
        value = term.is_zero
        if value is not None:
            return not bool(value)
    if term.is_Add and term.has(sp.Pow):
        sqrt_sum = square_root_sum_test(term)
        if sqrt_sum.verdict is Verdict.NONZERO_PROVEN:
            return True
        if sqrt_sum.verdict is Verdict.ZERO_PROVEN:
            return False
    if ops > cfg.EXP_INDEP_GAP_MAX_OPS:
        return None
    result = run_with_time_budget(
        algebraic_gap_test, term,
        seconds=cfg.EXP_INDEP_GAP_TIMEOUT, default=None,
    )
    if result is None:
        return None
    if result.verdict is Verdict.NONZERO_PROVEN:
        return True
    if result.verdict is Verdict.ZERO_PROVEN:
        return False
    return None


def _exp_term(term: sp.Expr):
    """Return ``(coefficient, exponent)`` for c*exp(a), constants as exp(0)."""
    if term is sp.E:
        return sp.Integer(1), sp.Integer(1)
    if term == -sp.E:
        return sp.Integer(-1), sp.Integer(1)
    coeff, rest = term.as_coeff_Mul(rational=False)
    if rest.func is sp.exp and len(rest.args) == 1:
        return coeff, rest.args[0]
    if term.func is sp.exp and len(term.args) == 1:
        return sp.Integer(1), term.args[0]
    # Algebraic constants are coefficients of exp(0).
    if _alg_coeff(term):
        return term, sp.Integer(0)
    if term.is_Mul:
        exp_parts = [arg for arg in term.args if arg.func is sp.exp and len(arg.args) == 1]
        if len(exp_parts) == 1:
            exp_part = exp_parts[0]
            coeff = sp.Mul(*(arg for arg in term.args if arg is not exp_part))
            return coeff, exp_part.args[0]
    return None


def exponential_independence_nonzero(term: sp.Expr) -> bool | None:
    """Prove a small algebraic linear combination of exponentials nonzero.

    Lindemann--Weierstrass implies that exponentials of distinct algebraic
    numbers are linearly independent over the algebraic numbers.  This helper
    deliberately recognizes only small, structurally explicit sums.
    """
    term = sp.sympify(term)
    parts = term.args if term.is_Add else (term,)
    if len(parts) < 2 or len(parts) > cfg.EXP_INDEP_MAX_TERMS:
        return None
    pairs = []
    saw_nonzero_coeff = False
    for part in parts:
        item = _exp_term(part)
        if item is None:
            return None
        coeff, exponent = item
        if not _alg_coeff(coeff) or number_kind(exponent) is not NumberKind.ALGEBRAIC:
            return None
        if coeff.is_zero is False or (coeff.is_Rational and coeff != 0):
            saw_nonzero_coeff = True
        elif coeff.is_zero is True:
            continue
        pairs.append((coeff, exponent))
    if not saw_nonzero_coeff or len(pairs) < 2:
        return None
    exponents = [exponent for _, exponent in pairs]
    for i, left in enumerate(exponents):
        for right in exponents[i + 1:]:
            if left == right:
                return None
            if _algebraic_nonzero(left - right) is not True:
                return None
    return True


def clear_exp_indep_cache() -> None:
    _algebraic_nonzero.cache_clear()
