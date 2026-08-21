from __future__ import annotations

from dataclasses import dataclass
from functools import cmp_to_key, total_ordering

import sympy as sp

from ..exact_arithmetic import compare_exact_reals
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
        try:
            roots = sorted(
                sp.real_roots(self.polynomial.as_expr()),
                key=cmp_to_key(compare_exact_reals),
            )
            return roots[self.root_index]
        except (NotImplementedError, sp.PolynomialError, ValueError):
            pass

        # SymPy cannot construct RootOf objects over every algebraic coefficient
        # domain.  For low-degree fibers, recover a native radical only when its
        # real identity is certified by this root's isolating interval.  This is
        # presentation-only: the interval remains the authoritative root identity.
        if self.polynomial.degree() <= 4:
            try:
                radicals = sp.roots(self.polynomial.as_expr(), self.variable)
                candidates = []
                for expr, multiplicity in radicals.items():
                    if expr.is_real is not True:
                        continue
                    left_cmp = compare_exact_reals(expr, self.interval.left)
                    right_cmp = compare_exact_reals(expr, self.interval.right)
                    if left_cmp >= 0 and right_cmp <= 0:
                        candidates.extend([sp.simplify(expr)] * int(multiplicity))
                if len(candidates) == 1:
                    return candidates[0]
            except (
                NotImplementedError,
                sp.PolynomialError,
                ValueError,
                TypeError,
            ):
                pass

        # Keep an exact ordered-root selector when no readable exact presentation
        # can be certified.
        from ..reconstruct.root_functions import root_of

        return root_of(
            self.polynomial.as_expr(),
            self.variable,
            sp.Integer(self.root_index),
        )

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
