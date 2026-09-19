"""Exact fast path for small rational linear combinations of square roots."""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from . import _config as cfg
from ._errors import EXACT_METHOD_ERRORS
from ._result import Verdict, ZeroClassification


@lru_cache(maxsize=cfg.NORMAL_FORM_CACHE_SIZE)
def _squarefree_int(value: int) -> tuple[sp.Integer, int] | None:
    """Return ``(outside, squarefree_part)`` for a small positive integer."""
    if value <= 0 or value.bit_length() > cfg.SQRT_SUM_RAD_BITS:
        return None
    try:
        factors = sp.factorint(value)
    except EXACT_METHOD_ERRORS:
        return None
    outside = 1
    squarefree = 1
    for prime, power in factors.items():
        outside *= int(prime) ** (int(power) // 2)
        if int(power) % 2:
            squarefree *= int(prime)
    return sp.Integer(outside), squarefree


def _sqrt_term(term: sp.Expr) -> tuple[int, sp.Rational] | None:
    """Return a squarefree radical key and rational coefficient."""
    term = sp.sympify(term)
    if term.is_Rational:
        return 1, sp.Rational(term)
    coeff, rest = term.as_coeff_Mul(rational=True)
    if not coeff.is_Rational:
        return None
    if not (rest.is_Pow and rest.exp == sp.Rational(1, 2)
            and rest.base.is_Integer and rest.base.is_positive):
        return None
    reduced = _squarefree_int(int(rest.base))
    if reduced is None:
        return None
    outside, key = reduced
    return key, sp.Rational(coeff) * outside


def square_root_sum_test(term: sp.Expr) -> ZeroClassification:
    """Decide small sums ``sum(q_i*sqrt(n_i))`` exactly.

    Distinct squarefree positive integer radicals are linearly independent
    over the rationals.  The method therefore reduces the expression to a
    sparse coefficient vector and never constructs a number field or minimal
    polynomial.
    """
    term = sp.sympify(term)
    if term.free_symbols or not term.is_Add:
        return ZeroClassification(Verdict.UNKNOWN, "sqrt-sum", detail="not a closed additive square-root sum")
    parts = term.args
    if len(parts) < 2 or len(parts) > cfg.SQRT_SUM_MAX_TERMS:
        return ZeroClassification(Verdict.UNKNOWN, "sqrt-sum", detail="square-root term budget exceeded")
    coeffs: dict[int, sp.Rational] = {}
    for part in parts:
        item = _sqrt_term(part)
        if item is None:
            return ZeroClassification(Verdict.UNKNOWN, "sqrt-sum", detail="unsupported square-root term")
        key, coeff = item
        coeffs[key] = coeffs.get(key, sp.Rational(0)) + coeff
    coeffs = {key: coeff for key, coeff in coeffs.items() if coeff != 0}
    if not coeffs:
        return ZeroClassification(
            Verdict.ZERO_PROVEN,
            "sqrt-sum",
            detail="squarefree radical coefficient vector cancels exactly",
            evidence="squarefree-radical-basis",
        )
    return ZeroClassification(
        Verdict.NONZERO_PROVEN,
        "sqrt-sum",
        detail="squarefree radical coefficient vector is nonzero",
        evidence="squarefree-radical-basis",
    )


def clear_sqrt_sum_cache() -> None:
    _squarefree_int.cache_clear()
