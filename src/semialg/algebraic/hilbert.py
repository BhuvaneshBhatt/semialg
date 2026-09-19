"""Exact Hilbert-series and degree invariants for polynomial ideals.

The implementation works from the leading monomial ideal of an exact,
degree-compatible Groebner basis.  Hilbert series are computed by the Taylor
(inclusion--exclusion) resolution of a monomial ideal, so the resulting
invariants are independent of any triangular decomposition search.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from itertools import combinations
from math import factorial

import sympy as sp

from .equality_ideal import EqualityIdealContext, _leading_exponents
from .rational_univariate.representation import RationalUnivariateError


@dataclass(frozen=True)
class HilbertData:
    """Exact Hilbert data of ``k[x]/I`` for a polynomial ideal ``I``."""

    variables: tuple[sp.Symbol, ...]
    dimension: int
    degree: int
    numerator: sp.Expr
    denominator_power: int
    leading_monomials: tuple[tuple[int, ...], ...]

    @property
    def series(self) -> sp.Expr:
        """Return the reduced Hilbert series in a fresh formal variable."""
        t = sp.Symbol("t")
        return sp.cancel(self.numerator / (1 - t) ** self.denominator_power)


def _divides(left: Sequence[int], right: Sequence[int]) -> bool:
    return all(int(a) <= int(b) for a, b in zip(left, right, strict=True))


def _minimal_monomial_generators(
    exponents: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    unique = sorted({tuple(map(int, exponent)) for exponent in exponents})
    return tuple(
        exponent
        for exponent in unique
        if not any(other != exponent and _divides(other, exponent) for other in unique)
    )


def _lcm_exponent(exponents: Sequence[Sequence[int]], variable_count: int) -> tuple[int, ...]:
    if not exponents:
        return (0,) * variable_count
    return tuple(max(int(exponent[i]) for exponent in exponents) for i in range(variable_count))


def _hilbert_numerator(
    monomials: Sequence[Sequence[int]], variable_count: int, t: sp.Symbol
) -> sp.Expr:
    """Return the Taylor-resolution numerator of a monomial ideal quotient."""
    generators = _minimal_monomial_generators(monomials)
    numerator = sp.Integer(0)
    for size in range(len(generators) + 1):
        sign = -1 if size % 2 else 1
        for subset in combinations(generators, size):
            exponent = _lcm_exponent(subset, variable_count)
            numerator += sign * t ** sum(exponent)
    return sp.Poly(sp.expand(numerator), t, domain=sp.ZZ).as_expr()


def _context(
    equations: Iterable[sp.Expr | sp.Equality], variables: Sequence[sp.Symbol]
) -> EqualityIdealContext:
    try:
        return EqualityIdealContext.build(equations, variables)
    except RationalUnivariateError as exc:
        raise NotImplementedError(
            "Hilbert invariants require exact polynomial coefficients"
        ) from exc


def ideal_hilbert_data(
    equations: Iterable[sp.Expr | sp.Equality], variables: Sequence[sp.Symbol]
) -> HilbertData:
    """Return exact Hilbert-series data and affine degree of ``<equations>``.

    A degree-compatible grevlex basis is used, so ``I`` and its initial
    monomial ideal have the same Hilbert function.  If
    ``H(t) = Q(t)/(1-t)^d`` is reduced at ``t=1``, the affine multiplicity
    (degree) is ``Q(1)``.
    """
    vars_ = tuple(variables)
    if not vars_:
        raise ValueError("Hilbert invariants require ambient variables")
    context = _context(equations, vars_)
    if context.inconsistent:
        return HilbertData(vars_, -1, 0, sp.Integer(0), 0, ((0,) * len(vars_),))
    if context.groebner_basis is None:
        leading: tuple[tuple[int, ...], ...] = tuple()
    else:
        leading = _minimal_monomial_generators(_leading_exponents(context.groebner_basis))
    t = sp.Symbol("t")
    numerator = _hilbert_numerator(leading, len(vars_), t)
    dimension = context.dimension
    cancellation = len(vars_) - dimension
    divisor = sp.Poly((1 - t) ** cancellation, t, domain=sp.ZZ)
    quotient, remainder = sp.div(sp.Poly(numerator, t, domain=sp.ZZ), divisor)
    if remainder.as_expr() != 0:
        raise ArithmeticError("Hilbert numerator did not have the expected dimension cancellation")
    reduced_numerator = sp.expand(quotient.as_expr())
    degree = int(sp.expand(reduced_numerator).subs(t, 1))
    if degree < 0:
        raise ArithmeticError("computed a negative ideal degree")
    return HilbertData(vars_, dimension, degree, reduced_numerator, dimension, leading)


def ideal_degree(equations: Iterable[sp.Expr | sp.Equality], variables: Sequence[sp.Symbol]) -> int:
    """Return the exact affine degree/multiplicity of a polynomial ideal."""
    return ideal_hilbert_data(equations, variables).degree


def hilbert_function(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
    degree: int,
) -> int:
    """Return ``dim_k (k[x]/in(I))_degree`` exactly."""
    if degree < 0:
        return 0
    data = ideal_hilbert_data(equations, variables)
    if data.dimension < 0:
        return 0
    generators = data.leading_monomials
    n = len(data.variables)
    value = 0
    for size in range(len(generators) + 1):
        sign = -1 if size % 2 else 1
        for subset in combinations(generators, size):
            shift = sum(_lcm_exponent(subset, n))
            residual = degree - shift
            if residual >= 0:
                value += sign * int(sp.binomial(residual + n - 1, n - 1))
    return int(value)


def hilbert_polynomial(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
    *,
    symbol: sp.Symbol | None = None,
) -> sp.Expr:
    """Return the eventual Hilbert polynomial of the graded quotient."""
    data = ideal_hilbert_data(equations, variables)
    d = sp.Symbol("d", integer=True, nonnegative=True) if symbol is None else symbol
    if data.dimension <= 0:
        return sp.Integer(0)
    n = len(data.variables)
    expression = sp.Integer(0)
    generators = data.leading_monomials
    for size in range(len(generators) + 1):
        sign = -1 if size % 2 else 1
        for subset in combinations(generators, size):
            shift = sum(_lcm_exponent(subset, n))
            term = sp.prod(d - shift + n - 1 - j for j in range(n - 1)) / factorial(n - 1)
            expression += sign * term
    return sp.expand(sp.cancel(expression))


__all__ = [
    "HilbertData",
    "hilbert_function",
    "hilbert_polynomial",
    "ideal_degree",
    "ideal_hilbert_data",
]
