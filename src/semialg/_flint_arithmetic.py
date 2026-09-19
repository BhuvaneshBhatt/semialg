"""Centralized exact FLINT/Arb acceleration for integer/rational polynomials.

This module is internal.  Public semialg objects remain SymPy/
semialg abstractions; FLINT values are short-lived arithmetic representations.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from math import gcd

import sympy as sp

try:
    from flint import arb, ctx, fmpq, fmpq_mpoly_ctx, fmpq_poly, fmpz_mpoly_ctx, fmpz_poly
except ImportError as exc:  # pragma: no cover - packaging guarantees the dependency
    raise ImportError("semialg requires python-flint>=0.9.0") from exc


@dataclass(frozen=True)
class FlintUnivariate:
    poly: object
    variable: sp.Symbol
    domain: str


@dataclass(frozen=True)
class FlintMultivariate:
    poly: object
    generators: tuple[sp.Symbol, ...]
    domain: str


def _as_fmpq(value: sp.Expr) -> fmpq:
    rational = sp.Rational(value)
    return fmpq(int(rational.p), int(rational.q))


def _poly_ground_domain(poly: sp.Poly) -> str | None:
    if poly.domain == sp.ZZ:
        return "ZZ"
    if poly.domain == sp.QQ:
        return "QQ"
    return None


def to_flint_univariate(poly: sp.Poly) -> FlintUnivariate | None:
    """Convert a univariate ZZ/QQ SymPy polynomial to a transient FLINT value."""

    if poly.is_multivariate:
        return None
    domain = _poly_ground_domain(poly)
    if domain is None:
        return None
    coeffs_high = poly.all_coeffs()
    coeffs_low = list(reversed(coeffs_high))
    if domain == "ZZ":
        value = fmpz_poly([int(c) for c in coeffs_low])
    else:
        value = fmpq_poly([_as_fmpq(c) for c in coeffs_low])
    return FlintUnivariate(value, poly.gens[0], domain)


def from_flint_univariate(value: object, variable: sp.Symbol, *, domain: str) -> sp.Poly:
    coeffs: list[sp.Expr] = []
    for coeff in value:  # type: ignore[operator]
        if domain == "ZZ":
            coeffs.append(sp.Integer(int(coeff)))
        else:
            coeffs.append(sp.Rational(int(coeff.p), int(coeff.q)))
    expr = sum(coeff * variable**index for index, coeff in enumerate(coeffs))
    return sp.Poly(expr, variable, domain=sp.ZZ if domain == "ZZ" else sp.QQ)


def to_flint_multivariate(poly: sp.Poly) -> FlintMultivariate | None:
    """Convert a multivariate ZZ/QQ polynomial without exposing FLINT upstream."""

    domain = _poly_ground_domain(poly)
    if domain is None:
        return None
    names = [f"x{index}" for index in range(len(poly.gens))]
    ctx_type = fmpz_mpoly_ctx if domain == "ZZ" else fmpq_mpoly_ctx
    flint_ctx = ctx_type.get(names)
    terms = {}
    for powers, coeff in poly.terms():
        terms[tuple(int(power) for power in powers)] = (
            int(coeff) if domain == "ZZ" else _as_fmpq(coeff)
        )
    return FlintMultivariate(flint_ctx.from_dict(terms), tuple(poly.gens), domain)


def from_flint_multivariate(
    value: object, generators: Sequence[sp.Symbol], *, domain: str
) -> sp.Poly:
    expr = sp.Integer(0)
    for powers, coeff in value.to_dict().items():  # type: ignore[attr-defined]
        if domain == "ZZ":
            sym_coeff: sp.Expr = sp.Integer(int(coeff))
        else:
            sym_coeff = sp.Rational(int(coeff.p), int(coeff.q))
        monomial = sym_coeff
        for generator, power in zip(generators, powers, strict=True):
            monomial *= generator ** int(power)
        expr += monomial
    return sp.Poly(expr, *generators, domain=sp.ZZ if domain == "ZZ" else sp.QQ)


def flint_gcd(poly: sp.Poly, other: sp.Poly) -> sp.Poly | None:
    left = to_flint_univariate(poly)
    right = to_flint_univariate(other)
    if left is None or right is None or left.domain != right.domain:
        return None
    result = left.poly.gcd(right.poly)  # type: ignore[attr-defined]
    return from_flint_univariate(result, left.variable, domain=left.domain)


def flint_eval(poly: sp.Poly, point: sp.Rational) -> sp.Rational | sp.Integer | None:
    converted = to_flint_univariate(poly)
    if converted is None:
        return None
    if converted.domain == "ZZ" and point.q == 1:
        return sp.Integer(int(converted.poly(int(point))))  # type: ignore[operator]
    value = converted.poly(_as_fmpq(point))  # type: ignore[operator]
    if hasattr(value, "p"):
        return sp.Rational(int(value.p), int(value.q))
    return sp.Integer(int(value))


def flint_rational_roots(poly: sp.Poly) -> tuple[tuple[sp.Rational, int], ...] | None:
    """Return exact rational roots and multiplicities for a ZZ/QQ polynomial."""

    converted = to_flint_univariate(poly)
    if converted is None:
        return None
    roots = []
    for value, multiplicity in converted.poly.roots():  # type: ignore[attr-defined]
        if hasattr(value, "p"):
            root = sp.Rational(int(value.p), int(value.q))
        else:
            root = sp.Rational(int(value))
        roots.append((root, int(multiplicity)))
    return tuple(roots)


def flint_resultant(left: sp.Poly, right: sp.Poly, variable: sp.Symbol) -> sp.Expr | None:
    """Compute a ZZ/QQ resultant in FLINT and convert the result back to SymPy."""

    generators = tuple(dict.fromkeys((*left.gens, *right.gens)))
    try:
        left_poly = sp.Poly(left.as_expr(), *generators)
        right_poly = sp.Poly(right.as_expr(), *generators)
    except (sp.PolynomialError, TypeError, ValueError):
        return None
    if left_poly.domain not in (sp.ZZ, sp.QQ) or right_poly.domain not in (sp.ZZ, sp.QQ):
        return None
    domain = "QQ" if sp.QQ in (left_poly.domain, right_poly.domain) else "ZZ"
    if domain == "QQ":
        left_poly = sp.Poly(left_poly.as_expr(), *generators, domain=sp.QQ)
        right_poly = sp.Poly(right_poly.as_expr(), *generators, domain=sp.QQ)
    converted_left = to_flint_multivariate(left_poly)
    converted_right = to_flint_multivariate(right_poly)
    if converted_left is None or converted_right is None:
        return None
    try:
        variable_index = generators.index(variable)
    except ValueError:
        return None
    result = converted_left.poly.resultant(converted_right.poly, variable_index)  # type: ignore[attr-defined]
    converted = from_flint_multivariate(result, generators, domain=domain)
    return sp.expand(converted.as_expr())


def flint_squarefree_part(poly: sp.Poly) -> sp.Poly | None:
    """Return the primitive squarefree part of a ZZ/QQ polynomial via FLINT."""

    converted = to_flint_multivariate(poly)
    if converted is None:
        return None
    try:
        _unit, factors = converted.poly.factor_squarefree()  # type: ignore[attr-defined]
    except (AttributeError, TypeError, ValueError, RuntimeError):
        return None
    if not factors:
        return poly.primitive()[1]
    result = factors[0][0]
    for factor, _multiplicity in factors[1:]:
        result *= factor
    sympy_poly = from_flint_multivariate(result, converted.generators, domain=converted.domain)
    primitive = sympy_poly.primitive()[1]
    if primitive.LC().could_extract_minus_sign():
        primitive = -primitive
    return primitive


@contextmanager
def arb_precision(bits: int) -> Iterator[None]:
    old = ctx.prec
    ctx.prec = max(int(bits), 32)
    try:
        yield
    finally:
        ctx.prec = old


def _arb_rational(value: sp.Rational) -> arb:
    return arb(int(value.p)) / int(value.q)


def arb_newton_proposal(
    poly: sp.Poly,
    left: sp.Rational,
    right: sp.Rational,
    *,
    precision_bits: Sequence[int] = (64, 128, 256),
    max_denominator_bits: int = 48,
) -> sp.Rational | None:
    """Use Arb interval arithmetic to propose a controlled-height Newton split.

    The returned rational is only a proposal.  Exact Descartes/Sturm predicates
    in the caller certify every accepted subdivision.
    """

    converted = to_flint_univariate(poly)
    if converted is None:
        return None
    derivative = converted.poly.derivative()  # type: ignore[attr-defined]
    midpoint = sp.Rational(left + right, 2)
    for bits in precision_bits:
        with arb_precision(bits):
            mid_ball = _arb_rational(midpoint)
            try:
                value = converted.poly(mid_ball)  # type: ignore[operator]
                slope = derivative(mid_ball)
            except (TypeError, ValueError, ArithmeticError):
                continue
            if not slope.is_finite() or slope.contains(0):
                continue
            proposal = mid_ball - value / slope
            if not proposal.is_finite():
                continue
            center = proposal.mid()
            # Decimal extraction is proposal-only; exact predicates certify it.
            try:
                approximate = sp.Rational(str(center))
            except (ValueError, TypeError):
                continue
            if not left < approximate < right:
                continue
            # Keep split heights controlled so exact transformed coefficients do
            # not explode.  limit_denominator is deterministic and proposal-only.
            numerator = int(approximate.p)
            denominator = int(approximate.q)
            if denominator.bit_length() > max_denominator_bits:
                scale = 1 << max_denominator_bits
                scaled_num = int(approximate.p) * scale
                scaled_den = int(approximate.q)
                quotient, remainder = divmod(abs(scaled_num), scaled_den)
                rounded = quotient + int(2 * remainder >= scaled_den)
                numerator = rounded if scaled_num >= 0 else -rounded
                common = gcd(abs(numerator), scale)
                approximate = sp.Rational(numerator // common, scale // common)
            if left < approximate < right:
                return approximate
    return None


__all__ = [
    "FlintMultivariate",
    "FlintUnivariate",
    "arb_newton_proposal",
    "flint_eval",
    "flint_gcd",
    "flint_rational_roots",
    "flint_resultant",
    "flint_squarefree_part",
    "from_flint_multivariate",
    "from_flint_univariate",
    "to_flint_multivariate",
    "to_flint_univariate",
]
