from __future__ import annotations

import sympy as sp

from .samples import AlgebraicRoot, RationalSample, Sample


def _as_sample(value: Sample | sp.Expr | None) -> Sample | None:
    if value is None:
        return None
    if isinstance(value, (RationalSample, AlgebraicRoot)):
        return value
    raise TypeError("CAD sector bounds must be explicit semialg Sample objects")


def choose_sector_sample(left: Sample | None, right: Sample | None) -> RationalSample:
    """Choose a rational sample in the sector between two explicit bounds."""

    left = _as_sample(left)
    right = _as_sample(right)
    if left is None and right is None:
        return RationalSample(sp.Integer(0))
    if left is None:
        assert right is not None
        return RationalSample(right.interval.left - 1)
    if right is None:
        return RationalSample(left.interval.right + 1)
    if left.interval.right < right.interval.left:
        return RationalSample(sp.Rational(left.interval.right + right.interval.left, 2))
    # Overlapping isolating intervals can occur before refinement.  Establish
    # the order exactly, then refine algebraic isolating intervals until a
    # rational separator is available.
    from .comparison import compare_samples
    from .roots import refine_isol_intv

    if compare_samples(left, right) >= 0:
        raise ValueError("left sector bound must be strictly smaller than right sector bound")
    refined_left, refined_right = left, right
    for _ in range(32):
        if refined_left.interval.right < refined_right.interval.left:
            return RationalSample(
                sp.Rational(refined_left.interval.right + refined_right.interval.left, 2)
            )
        if isinstance(refined_left, AlgebraicRoot):
            refined_left = refine_isol_intv(refined_left, steps=2)
        if isinstance(refined_right, AlgebraicRoot):
            refined_right = refine_isol_intv(refined_right, steps=2)
    raise ValueError("could not obtain disjoint exact isolating intervals for sector bounds")
