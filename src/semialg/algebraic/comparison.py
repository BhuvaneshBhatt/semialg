from __future__ import annotations

import sympy as sp

from .cache import CACHE, sample_expr_key
from .roots import refine_isol_intv
from .samples import AlgebraicRoot, RationalSample, Sample, sample_to_expr


def compare_samples(left: Sample, right: Sample) -> int:
    """Compare two explicit sample objects exactly enough for CAD ordering."""

    if not isinstance(left, (RationalSample, AlgebraicRoot)) or not isinstance(
        right, (RationalSample, AlgebraicRoot)
    ):
        raise TypeError("compare_samples requires explicit semialg Sample objects")
    left_key = sample_expr_key(left)
    right_key = sample_expr_key(right)
    canonical = (left_key, right_key) if left_key <= right_key else (right_key, left_key)
    orientation = 1 if left_key <= right_key else -1
    cached = CACHE.comparisons.get(canonical)
    if cached is not None:
        CACHE.stats.comparison_hits += 1
        return orientation * cached
    CACHE.stats.comparison_misses += 1
    if isinstance(left, RationalSample) and isinstance(right, RationalSample):
        result = int(bool(left.value > right.value)) - int(bool(left.value < right.value))
        CACHE.comparisons.put(canonical, orientation * result)
        return result

    interval_order = left.interval.strict_order(right.interval)
    if interval_order is not None:
        CACHE.comparisons.put(canonical, orientation * interval_order)
        return interval_order
    if left.interval.right == right.interval.left:
        point = left.interval.right
        left_at_point = (isinstance(left, RationalSample) and left.value == point) or (
            isinstance(left, AlgebraicRoot) and sp.simplify(left.polynomial.eval(point)) == 0
        )
        right_at_point = (isinstance(right, RationalSample) and right.value == point) or (
            isinstance(right, AlgebraicRoot) and sp.simplify(right.polynomial.eval(point)) == 0
        )
        result = 0 if left_at_point and right_at_point else -1
        CACHE.comparisons.put(canonical, orientation * result)
        return result
    if right.interval.right == left.interval.left:
        point = right.interval.right
        left_at_point = (isinstance(left, RationalSample) and left.value == point) or (
            isinstance(left, AlgebraicRoot) and sp.simplify(left.polynomial.eval(point)) == 0
        )
        right_at_point = (isinstance(right, RationalSample) and right.value == point) or (
            isinstance(right, AlgebraicRoot) and sp.simplify(right.polynomial.eval(point)) == 0
        )
        result = 0 if left_at_point and right_at_point else 1
        CACHE.comparisons.put(canonical, orientation * result)
        return result

    diff = sp.simplify(sample_to_expr(left) - sample_to_expr(right))
    if diff == 0:
        CACHE.comparisons.put(canonical, 0)
        return 0
    sign = sp.sign(diff)
    if sign in (-1, 1):
        result = int(sign)
        CACHE.comparisons.put(canonical, orientation * result)
        return result
    # Last resort: refine disjoint intervals before using high-precision sign.
    refined_left = left
    refined_right = right
    for _ in range(64):
        if isinstance(refined_left, AlgebraicRoot):
            refined_left = refine_isol_intv(refined_left, steps=2)
        if isinstance(refined_right, AlgebraicRoot):
            refined_right = refine_isol_intv(refined_right, steps=2)
        interval_order = refined_left.interval.strict_order(refined_right.interval)
        if interval_order is not None:
            CACHE.comparisons.put(canonical, orientation * interval_order)
            return interval_order
        CACHE.stats.comparison_refinements += 1
    if isinstance(refined_left, AlgebraicRoot) and isinstance(refined_right, AlgebraicRoot):
        if refined_left.root_expr is None or refined_right.root_expr is None:
            try:
                var = refined_left.variable
                right_expr = refined_right.polynomial.as_expr().subs(refined_right.variable, var)
                left_poly = sp.Poly(refined_left.polynomial.as_expr(), var, extension=True)
                right_poly = sp.Poly(right_expr, var, extension=True)
                common = sp.gcd(left_poly, right_poly)
                overlap_left = max(refined_left.interval.left, refined_right.interval.left)
                overlap_right = min(refined_left.interval.right, refined_right.interval.right)
                if (
                    common.degree() > 0
                    and overlap_left <= overlap_right
                    and int(common.count_roots(overlap_left, overlap_right)) > 0
                ):
                    CACHE.comparisons.put(canonical, 0)
                    return 0
            except (
                NotImplementedError,
                sp.PolynomialError,
                sp.polys.polyerrors.DomainError,
                ValueError,
                TypeError,
            ):
                pass

    # Do not turn a fixed-precision approximation into an "exact" CAD
    # ordering.  Try a canonical algebraic-number representation; otherwise
    # fail conservatively so callers can use another exact strategy.
    try:
        algebraic = sp.polys.numberfields.to_number_field(diff).to_root()
        if algebraic == 0 or algebraic.is_zero is True:
            result = 0
        elif algebraic.is_positive is True:
            result = 1
        elif algebraic.is_negative is True:
            result = -1
        else:
            exact_sign = sp.sign(algebraic)
            if exact_sign not in (-1, 0, 1):
                raise ValueError("exact algebraic comparison remained undecidable")
            result = int(exact_sign)
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError) as exc:
        raise ValueError("could not compare algebraic samples exactly") from exc
    CACHE.comparisons.put(canonical, orientation * result)
    return result


def sort_samples(samples: list[Sample] | tuple[Sample, ...]) -> tuple[Sample, ...]:
    out = list(samples)
    for i in range(1, len(out)):
        item = out[i]
        pos = i
        while pos > 0 and compare_samples(out[pos - 1], item) > 0:
            out[pos] = out[pos - 1]
            pos -= 1
        out[pos] = item
    unique: list[Sample] = []
    for sample in out:
        if not unique or compare_samples(unique[-1], sample) != 0:
            unique.append(sample)
    return tuple(unique)
