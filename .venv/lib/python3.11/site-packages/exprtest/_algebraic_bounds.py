"""Certified separation bounds for exact algebraic values."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import sympy as sp

from . import _config as cfg
from ._algebraic_model import (
    _model_denominator_is_proven_nonzero,
    build_algebraic_model,
    classify_algebraic_expression,
)
from ._cost import stage_allowed, within_budget
from ._errors import EXACT_METHOD_ERRORS
from ._memo import minpoly_for
from ._result import Verdict, ZeroClassification


@dataclass(frozen=True)
class AlgebraicGap:
    """Certified lower bound for the magnitude of a nonzero algebraic value."""

    lower: sp.Expr
    polynomial: sp.Poly | None
    method: str = "reciprocal-cauchy"


def _ceil_root(value: sp.Rational, degree: int) -> sp.Integer:
    """Smallest integer n with n**degree >= a nonnegative rational value."""
    value = sp.Rational(value)
    if value <= 1:
        return sp.Integer(1)
    low = 1
    high = 2
    while sp.Integer(high) ** degree < value:
        high *= 2
    while low + 1 < high:
        mid = (low + high) // 2
        if sp.Integer(mid) ** degree >= value:
            high = mid
        else:
            low = mid
    return sp.Integer(high)


def _cauchy_gap(poly: sp.Poly) -> sp.Rational | None:
    coeffs = poly.all_coeffs()
    constant = coeffs[-1]
    if constant == 0:
        return None
    ratio = max((abs(sp.Rational(value, constant)) for value in coeffs[:-1]),
                default=sp.Integer(0))
    return sp.Rational(1, 1) / (1 + ratio)


def _fujiwara_gap(poly: sp.Poly) -> sp.Rational | None:
    """Rationalized reciprocal Fujiwara lower root bound."""
    coeffs = poly.all_coeffs()
    constant = coeffs[-1]
    degree = poly.degree()
    if constant == 0 or degree < 1:
        return None
    best = sp.Integer(1)
    for offset in range(1, degree + 1):
        coeff = coeffs[-1 - offset]
        if coeff == 0:
            continue
        ratio = abs(sp.Rational(coeff, constant))
        best = max(best, _ceil_root(ratio, offset))
    return sp.Rational(1, 2 * int(best))


def _dominant_radius(poly: sp.Poly, radius: sp.Rational) -> bool:
    """Whether the leading term strictly dominates all others at ``radius``."""
    if radius <= 0:
        return False
    coeffs = poly.all_coeffs()
    degree = poly.degree()
    lead = abs(sp.Rational(coeffs[0])) * radius**degree
    tail = sp.Rational(0)
    for offset, coeff in enumerate(coeffs[1:], 1):
        power = degree - offset
        tail += abs(sp.Rational(coeff)) * radius**power
    return lead > tail


def _dyadic_root_upper(poly: sp.Poly) -> sp.Rational:
    """Certified root-radius bound refined by exact dyadic bisection.

    Once the leading term dominates the coefficient tail on ``|z| = R``,
    Rouché's theorem places every root strictly inside that circle. Starting
    from Cauchy/Fujiwara, a few exact bisection steps often tighten the radius
    substantially without root isolation or floating-point arithmetic.
    """
    high = _root_upper(poly)
    if high <= 0:
        high = sp.Integer(1)
    while not _dominant_radius(poly, high):
        high *= 2
    low = sp.Rational(0)
    for _ in range(cfg.ALG_GAP_DYADIC_STEPS):
        mid = (low + high) / 2
        if _dominant_radius(poly, mid):
            high = mid
        else:
            low = mid
    return sp.Rational(high)


def _dyadic_gap(poly: sp.Poly) -> sp.Rational | None:
    """Lower bound nonzero roots via a refined reciprocal root radius."""
    coeffs = poly.all_coeffs()
    if not coeffs or coeffs[-1] == 0:
        return None
    var = poly.gens[0]
    reciprocal = sp.Poly.from_list(list(reversed(coeffs)), gens=var, domain=sp.QQ)
    upper = _dyadic_root_upper(reciprocal)
    return None if upper <= 0 else sp.Rational(1, 1) / upper


@lru_cache(maxsize=cfg.EXACT_BOUND_CACHE_SIZE)
def _poly_gap_cheap_cached(term: sp.Expr) -> AlgebraicGap | None:
    """Return the best inexpensive minimal-polynomial separation bound."""
    marker = sp.Dummy("z")
    try:
        poly = minpoly_for(term, marker).clear_denoms()[1]
        bounds = [("reciprocal-cauchy", _cauchy_gap(poly)),
                  ("reciprocal-fujiwara", _fujiwara_gap(poly))]
        valid = [(name, bound) for name, bound in bounds
                 if bound is not None and bound > 0]
        if not valid:
            return None
        name, lower = max(valid, key=lambda item: item[1])
        return AlgebraicGap(lower, poly, name)
    except EXACT_METHOD_ERRORS:
        return None


@lru_cache(maxsize=cfg.EXACT_BOUND_CACHE_SIZE)
def _poly_gap_cached(term: sp.Expr) -> AlgebraicGap | None:
    """Refine a cheap bound only when a caller explicitly asks for strength."""
    cheap = _poly_gap_cheap_cached(term)
    if cheap is None:
        return None
    try:
        dyadic = _dyadic_gap(cheap.polynomial) if cheap.polynomial is not None else None
        if dyadic is not None and dyadic > cheap.lower:
            return AlgebraicGap(dyadic, cheap.polynomial, "reciprocal-dyadic")
    except EXACT_METHOD_ERRORS:
        pass
    return cheap


def _root_upper(poly: sp.Poly) -> sp.Rational:
    """Sharper rational upper bound for every root of ``poly``."""
    coeffs = poly.all_coeffs()
    lead = coeffs[0]
    tail = [abs(sp.Rational(value, lead)) for value in coeffs[1:]]
    cauchy = sp.Rational(1) + max(tail, default=sp.Integer(0))
    degree = poly.degree()
    best = sp.Integer(0)
    for offset in range(1, degree + 1):
        coeff = coeffs[offset]
        if coeff == 0:
            continue
        ratio = abs(sp.Rational(coeff, lead))
        best = max(best, _ceil_root(ratio, offset))
    fujiwara = 2 * best if best else sp.Integer(1)
    return min(cauchy, sp.Rational(fujiwara))


def _tree_upper(term: sp.Expr, limits: dict[sp.Symbol, sp.Rational]) -> sp.Rational | None:
    """Bound an expression tree over algebraic generators without expansion."""
    if term.is_Rational:
        return abs(sp.Rational(term))
    if term.is_Symbol and term in limits:
        return limits[term]
    if term.is_Add:
        parts = [_tree_upper(arg, limits) for arg in term.args]
        if any(part is None for part in parts):
            return None
        return sum((part for part in parts if part is not None), sp.Rational(0))
    if term.is_Mul:
        value = sp.Rational(1)
        for arg in term.args:
            part = _tree_upper(arg, limits)
            if part is None:
                return None
            value *= part
        return value
    if term.is_Pow and term.exp.is_Integer and term.exp.is_nonnegative:
        base = _tree_upper(term.base, limits)
        return None if base is None else base ** int(term.exp)
    return None


def _poly_upper(term: sp.Expr, vars_: tuple[sp.Symbol, ...],
                limits: tuple[sp.Rational, ...]) -> sp.Rational:
    """Exact L1 upper bound for a rational-coefficient polynomial."""
    poly = sp.Poly(term, *vars_, domain=sp.QQ)
    total = sp.Rational(0)
    for powers, coeff in poly.terms():
        size = abs(sp.Rational(coeff))
        for limit, power in zip(limits, powers):
            size *= limit ** power
        total += size
    return total


@lru_cache(maxsize=cfg.EXACT_BOUND_CACHE_SIZE)
def _resultant_gap_cached(term: sp.Expr, include_poly: bool = True) -> AlgebraicGap | None:
    """Bound one conjugate using an iterated exact resultant product.

    For monic defining polynomials the iterated resultant is the product of
    the numerator over every tuple of generator conjugates. A nonzero product
    has exact rational magnitude, while all other factors are bounded above by
    an L1 polynomial bound. This yields a certified lower bound for the
    selected algebraic value.
    """
    model = build_algebraic_model(term)
    if model is None or not model.generators:
        return None
    degree_product = model.degree_product
    if degree_product > cfg.ALG_RESULTANT_MAX_DEGREE or not stage_allowed(term, "resultant"):
        return None
    if sp.count_ops(model.numerator) > cfg.ALG_RESULTANT_MAX_OPS:
        return None
    try:
        monic = tuple(poly.monic() for poly in model.minpolys)
        limits = tuple(_dyadic_root_upper(poly) for poly in monic)
        limit_map = dict(zip(model.variables, limits))
        num_poly = _poly_upper(model.numerator, model.variables, limits)
        den_poly = _poly_upper(model.denominator, model.variables, limits)
        num_tree = _tree_upper(model.numerator, limit_map)
        den_tree = _tree_upper(model.denominator, limit_map)
        num_upper = min(num_poly, num_tree) if num_tree is not None else num_poly
        den_upper = min(den_poly, den_tree) if den_tree is not None else den_poly
        if num_upper <= 0 or den_upper <= 0:
            return None
        product = model.numerator
        for var, rel in zip(model.variables, monic):
            product = sp.resultant(rel.as_expr(), product, var)
            if sp.count_ops(product) > cfg.ALG_RESULTANT_MAX_OPS:
                return None
        product = sp.cancel(product)
        if product.free_symbols or product == 0 or not product.is_Rational:
            return None
        if not _model_denominator_is_proven_nonzero(model):
            return None
        lower_num = abs(sp.Rational(product)) / (num_upper ** (degree_product - 1))
        lower = sp.cancel(lower_num / den_upper)
        if lower <= 0:
            return None
        poly = None
        if include_poly:
            marker = sp.Dummy("z")
            poly = minpoly_for(term, marker).clear_denoms()[1]
        return AlgebraicGap(lower, poly, "conjugate-resultant")
    except EXACT_METHOD_ERRORS:
        return None


def algebraic_gap_bound(term: sp.Expr, refine: bool = True) -> AlgebraicGap | None:
    """Return a certified algebraic separation bound.

    ``refine=False`` stops after the inexpensive Cauchy/Fujiwara layer. The
    fast zero oracle uses that layer first because *any* positive certified
    lower bound is already enough to prove nonzero. Stronger dyadic and
    resultant bounds are computed only for callers that actually need them.
    """
    if not within_budget(sp.sympify(term)):
        return None

    term = sp.sympify(term)
    info = classify_algebraic_expression(term)
    if not info.is_algebraic:
        return None
    cheap = _poly_gap_cheap_cached(term)
    if cheap is None or not refine:
        return cheap
    bounds = [bound for bound in (_poly_gap_cached(term), _resultant_gap_cached(term))
              if bound is not None]
    return max(bounds, key=lambda item: item.lower) if bounds else cheap


def algebraic_gap_test(term: sp.Expr) -> ZeroClassification:
    """Use exact minimal-polynomial and separation bounds to decide zero."""
    if not within_budget(sp.sympify(term)):
        return ZeroClassification(Verdict.UNKNOWN, "algebraic-gap", detail="exact-method budget exceeded")

    term = sp.sympify(term)
    info = classify_algebraic_expression(term)
    if not info.is_algebraic:
        return ZeroClassification(Verdict.UNKNOWN, "algebraic-gap", detail=info.detail)
    marker = sp.Dummy("z")
    try:
        min_poly = minpoly_for(term, marker).clear_denoms()[1]
        if min_poly.degree() == 1 and min_poly.eval(0) == 0:
            return ZeroClassification(
                Verdict.ZERO_PROVEN, "algebraic-gap",
                detail="minimal polynomial is linear with root zero",
                evidence="minimal-polynomial",
            )
        # Any positive certified bound proves nonzero. Avoid the stronger
        # dyadic/resultant work unless a caller explicitly requests it.
        gap = algebraic_gap_bound(term, refine=False)
        if gap is not None:
            return ZeroClassification(
                Verdict.NONZERO_PROVEN, "algebraic-gap",
                detail=f"{gap.method} bound gives |value| >= {gap.lower}",
                evidence="algebraic-separation",
            )
    except EXACT_METHOD_ERRORS as exc:
        return ZeroClassification(Verdict.UNKNOWN, "algebraic-gap", detail=str(exc))
    return ZeroClassification(Verdict.UNKNOWN, "algebraic-gap", detail="no separation bound obtained")




def clear_bound_cache() -> None:
    """Clear cached minimal-polynomial and resultant separation bounds."""
    _poly_gap_cheap_cached.cache_clear()
    _poly_gap_cached.cache_clear()
    _resultant_gap_cached.cache_clear()
