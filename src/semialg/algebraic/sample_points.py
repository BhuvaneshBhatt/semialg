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
    # Overlapping isolating intervals should be rare after root sorting/deduping,
    # but choose a certified rational using exact expressions when necessary.
    left_value = sp.N(left.as_expr(), 100)
    right_value = sp.N(right.as_expr(), 100)
    if not left_value < right_value:
        raise ValueError("left sector bound must be strictly smaller than right sector bound")
    return RationalSample(sp.Rational(str((left_value + right_value) / 2)))
