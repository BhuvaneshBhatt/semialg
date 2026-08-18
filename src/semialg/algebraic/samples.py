from __future__ import annotations

from dataclasses import dataclass
from functools import total_ordering

import sympy as sp

from .intervals import RationalInterval


@total_ordering
@dataclass(frozen=True)
class RationalSample:
    """Exact rational CAD sample point."""

    value: sp.Rational

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", sp.Rational(self.value))

    @property
    def interval(self) -> RationalInterval:
        return RationalInterval(self.value, self.value)

    def as_expr(self) -> sp.Expr:
        return self.value

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, (RationalSample, AlgebraicRoot)):
            return NotImplemented
        from .comparison import compare_samples

        return compare_samples(self, other) < 0

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, (RationalSample, AlgebraicRoot)):
            return False
        from .comparison import compare_samples

        return compare_samples(self, other) == 0


@total_ordering
@dataclass(frozen=True)
class AlgebraicRoot:
    """A real algebraic root represented by a polynomial, expression, and isolating interval."""

    polynomial: sp.Poly
    interval: RationalInterval
    root_index: int
    multiplicity: int = 1
    root_expr: sp.Expr | None = None

    def __post_init__(self) -> None:
        if self.polynomial.is_multivariate:
            raise ValueError("AlgebraicRoot requires a univariate polynomial")
        if self.polynomial.degree() <= 0:
            raise ValueError("AlgebraicRoot requires a positive-degree polynomial")
        if self.multiplicity < 1:
            raise ValueError("multiplicity must be positive")

    @property
    def variable(self) -> sp.Symbol:
        return self.polynomial.gens[0]

    def as_expr(self) -> sp.Expr:
        if self.root_expr is not None:
            return self.root_expr
        roots = sorted(sp.real_roots(self.polynomial.as_expr()), key=lambda root: sp.N(root, 100))
        return roots[self.root_index]

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, (RationalSample, AlgebraicRoot)):
            return NotImplemented
        from .comparison import compare_samples

        return compare_samples(self, other) < 0

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, (RationalSample, AlgebraicRoot)):
            return False
        from .comparison import compare_samples

        return compare_samples(self, other) == 0


Sample = RationalSample | AlgebraicRoot


def sample_to_expr(sample: Sample) -> sp.Expr:
    if not isinstance(sample, (RationalSample, AlgebraicRoot)):
        raise TypeError("expected a semialg algebraic Sample object")
    return sample.as_expr()
