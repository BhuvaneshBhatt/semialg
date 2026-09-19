from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True, order=True)
class RationalInterval:
    """Closed rational interval used to isolate one real algebraic number."""

    left: sp.Rational
    right: sp.Rational

    def __post_init__(self) -> None:
        left = sp.Rational(self.left)
        right = sp.Rational(self.right)
        if left > right:
            raise ValueError("interval left endpoint must not exceed right endpoint")
        object.__setattr__(self, "left", left)
        object.__setattr__(self, "right", right)

    @property
    def width(self) -> sp.Rational:
        return self.right - self.left

    @property
    def midpoint(self) -> sp.Rational:
        return sp.Rational(self.left + self.right, 2)

    def contains_rational(self, value: sp.Rational) -> bool:
        value = sp.Rational(value)
        return self.left <= value <= self.right

    def is_point(self) -> bool:
        return self.left == self.right

    def is_disjoint_from(self, other: RationalInterval) -> bool:
        return self.right < other.left or other.right < self.left

    def strict_order(self, other: RationalInterval) -> int | None:
        """Return -1/1 when intervals are disjoint and ordered, else None."""

        if self.right < other.left:
            return -1
        if other.right < self.left:
            return 1
        return None

    def expand(self, radius: sp.Rational) -> RationalInterval:
        radius = abs(sp.Rational(radius))
        return RationalInterval(self.left - radius, self.right + radius)
