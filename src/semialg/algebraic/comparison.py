from __future__ import annotations

import sympy as sp

from .cache import CACHE
from .roots import refine_isol_intv
from .samples import AlgebraicRoot, RationalSample, Sample, sample_to_expr


def compare_samples(left: Sample, right: Sample) -> int:
    """Compare two explicit sample objects exactly enough for CAD ordering."""

    if not isinstance(left, (RationalSample, AlgebraicRoot)) or not isinstance(
        right, (RationalSample, AlgebraicRoot)
    ):
        raise TypeError("compare_samples requires explicit semialg Sample objects")
    if isinstance(left, RationalSample) and isinstance(right, RationalSample):
        return int(bool(left.value > right.value)) - int(bool(left.value < right.value))

    interval_order = left.interval.strict_order(right.interval)
    if interval_order is not None:
        return interval_order

    diff = sp.simplify(sample_to_expr(left) - sample_to_expr(right))
    if diff == 0:
        return 0
    sign = sp.sign(diff)
    if sign in (-1, 1):
        return int(sign)
    # Last resort: refine disjoint intervals before using high-precision sign.
    refined_left = left
    refined_right = right
    for _ in range(8):
        if isinstance(refined_left, AlgebraicRoot):
            refined_left = refine_isol_intv(refined_left, steps=2)
        if isinstance(refined_right, AlgebraicRoot):
            refined_right = refine_isol_intv(refined_right, steps=2)
        interval_order = refined_left.interval.strict_order(refined_right.interval)
        if interval_order is not None:
            return interval_order
        CACHE.stats.comparison_refinements += 1
    value = sp.N(diff, 160)
    return 1 if value > 0 else -1 if value < 0 else 0


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
