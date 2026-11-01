"""Canonical polynomial helpers shared by CAD projection and lifting."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import lru_cache
from typing import TypeAlias

import sympy as sp

from .._flint_arithmetic import flint_resultant
from ..errors import ExactComputationUnavailable

PolynomialKey: TypeAlias = tuple[
    tuple[sp.Symbol, ...], object, tuple[tuple[tuple[int, ...], sp.Expr], ...]
]


def _cad_input_polynomials(polynomials, variables):
    """Bind a real polynomial family to its full ambient coordinate tuple.

    A missing coordinate must not turn an input polynomial into a coefficient
    that normalization discards as a unit. Constants in zero coordinates are
    retained as scalar expressions, since SymPy Poly needs a generator.
    """
    variables = tuple(variables)
    if any(not isinstance(variable, sp.Symbol) for variable in variables):
        raise TypeError("CAD variables must be symbols")
    if len(set(variables)) != len(variables):
        raise ValueError("CAD variables must be distinct")
    ambient = set(variables)
    result = []
    for item in polynomials:
        expression = item.as_expr() if isinstance(item, sp.Poly) else sp.sympify(item)
        missing = expression.free_symbols - ambient
        if missing:
            raise ValueError(
                "CAD polynomial contains unassigned variables: "
                + ", ".join(str(v) for v in sorted(missing, key=sp.default_sort_key))
            )
        if not variables:
            from ..exact_arithmetic import exact_sign

            if expression.is_real is False:
                raise ValueError("CAD requires real polynomial coefficients")
            exact_sign(expression)
            result.append(expression)
            continue
        polynomial = (
            item
            if isinstance(item, sp.Poly) and item.gens == variables
            else sp.Poly(expression, *variables)
        )
        if any(coefficient.is_real is False for coefficient in polynomial.coeffs()):
            raise ValueError("CAD requires real polynomial coefficients")
        result.append(polynomial)
    return tuple(result)


def polynomial_structural_key(poly: sp.Poly) -> PolynomialKey:
    """Return a hashable exact identity key without expression serialization."""

    return (tuple(poly.gens), poly.domain, tuple(poly.terms()))


def polynomial_sort_key(poly: sp.Poly) -> tuple[object, tuple[object, ...]]:
    """Return a deterministic structural ordering key for normalized polynomials."""

    return (
        sp.default_sort_key(poly.as_expr()),
        tuple(sp.default_sort_key(gen) for gen in poly.gens),
    )


def polynomial_family_key(polys: Iterable[sp.Poly]) -> tuple[PolynomialKey, ...]:
    """Return an order-independent structural key for a polynomial family."""

    ordered = sorted(polys, key=polynomial_sort_key)
    return tuple(polynomial_structural_key(poly) for poly in ordered)


@lru_cache(maxsize=4096)
def polynomial_key(poly: sp.Poly) -> str:
    """Return an identity-preserving canonical expression key for CAD metadata.

    ``srepr`` retains Symbol assumptions, unlike pretty-printer text, so two
    same-name generators with different identities cannot collide in sign or
    provenance maps. The result remains cached because projection polynomials
    recur across many cells during lifting.
    """

    expr = poly.as_expr()
    readable = sp.sstr(expr)
    if all(sp.srepr(symbol) == sp.srepr(sp.Symbol(symbol.name)) for symbol in expr.free_symbols):
        return readable
    return readable + "\x1f" + sp.srepr(expr)


@dataclass(frozen=True)
class CoefficientDomainPolicy:
    """Cost guard for exact algebraic coefficient-domain construction.

    Natural exact ground domains (notably ``ZZ`` and ``QQ``) are always kept.
    Explicit algebraic number fields are constructed only when the extension
    degree and coefficient-expression complexity are small enough to plausibly
    beat SymPy's expression domain.  This is a performance policy only: every
    selected representation is exact and callers may always fall back to EX.
    """

    max_extension_degree: int = 4
    max_algebraic_generators: int = 2
    max_coefficient_ops: int = 80


@dataclass(frozen=True)
class CoefficientDomainDecision:
    domain: object
    strategy: str
    extension_degree: int | None = None
    algebraic_generators: int = 0
    coefficient_ops: int = 0


DEFAULT_DOMAIN_POLICY = CoefficientDomainPolicy()


def _algebraic_generator_count(expr: sp.Expr) -> int:
    generators: set[sp.Expr] = set()
    for atom in sp.preorder_traversal(expr):
        if (
            isinstance(atom, sp.Pow)
            and atom.exp.is_Rational
            and not atom.exp.is_Integer
            or getattr(atom, "is_AlgebraicNumber", False)
        ):
            generators.add(atom)
    return len(generators)


def _univariate_domain_data(
    expr: sp.Expr,
    var: sp.Symbol,
    *,
    algebraic_extension: bool = True,
    policy: CoefficientDomainPolicy = DEFAULT_DOMAIN_POLICY,
) -> tuple[CoefficientDomainDecision, sp.Poly | None]:
    """Select the domain and retain the polynomial that established the choice."""

    value = sp.sympify(expr)
    try:
        natural = sp.Poly(value, var)
    except (sp.PolynomialError, TypeError, ValueError):
        natural = None
    if (
        natural is not None
        and natural.domain != sp.EX
        and getattr(natural.domain, "is_Exact", False)
    ):
        return CoefficientDomainDecision(natural.domain, "natural_exact"), natural
    ops = int(sp.count_ops(value))
    generators = _algebraic_generator_count(value)
    if (
        not algebraic_extension
        or generators > policy.max_algebraic_generators
        or ops > policy.max_coefficient_ops
    ):
        return CoefficientDomainDecision(sp.EX, "expression_domain", None, generators, ops), None
    try:
        extended = sp.Poly(value, var, extension=True)
    except (sp.PolynomialError, TypeError, ValueError):
        return CoefficientDomainDecision(sp.EX, "expression_domain", None, generators, ops), None
    if extended.domain == sp.EX or not getattr(extended.domain, "is_Exact", False):
        return CoefficientDomainDecision(sp.EX, "expression_domain", None, generators, ops), None
    degree = None
    try:
        degree = int(extended.domain.ext.minpoly.degree())
    except (AttributeError, TypeError, ValueError):
        degree = 1
    if degree is not None and degree > policy.max_extension_degree:
        return CoefficientDomainDecision(
            sp.EX, "extension_too_large", degree, generators, ops
        ), None
    return (
        CoefficientDomainDecision(
            extended.domain, "small_algebraic_extension", degree, generators, ops
        ),
        extended,
    )


def coefficient_domain_decision(
    expr: sp.Expr,
    var: sp.Symbol,
    *,
    algebraic_extension: bool = True,
    policy: CoefficientDomainPolicy = DEFAULT_DOMAIN_POLICY,
) -> CoefficientDomainDecision:
    """Choose an exact univariate coefficient representation by cost policy."""
    return _univariate_domain_data(
        expr, var, algebraic_extension=algebraic_extension, policy=policy
    )[0]


def exact_univariate_poly(
    expr: sp.Expr,
    var: sp.Symbol,
    *,
    algebraic_extension: bool = True,
    policy: CoefficientDomainPolicy = DEFAULT_DOMAIN_POLICY,
) -> sp.Poly:
    """Return an exact univariate polynomial using selective domain promotion."""

    value = sp.sympify(expr)
    _decision, poly = _univariate_domain_data(
        value, var, algebraic_extension=algebraic_extension, policy=policy
    )
    return poly if poly is not None else sp.Poly(value, var, domain="EX")


@lru_cache(maxsize=1024)
def normalize_poly(poly: sp.Poly) -> sp.Poly | None:
    """Normalize a nonconstant polynomial up to content, sign, and multiplicity.

    Projection only needs the product of the distinct irreducible factors.  That
    is exactly the primitive square-free part, so avoid a full factorization here.
    """

    if poly.is_zero:
        return None
    primitive = poly.primitive()[1]
    if primitive.is_zero or primitive.total_degree() <= 0:
        return None
    try:
        result = primitive.sqf_part()
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        # Some exact extension/expression domains do not implement multivariate
        # square-free decomposition.  Preserve the previous factor-based fallback
        # for those domains rather than weakening normalization semantics.
        try:
            factors = sp.factor_list(primitive.as_expr(), *primitive.gens)[1]
        except (sp.PolynomialError, ValueError, TypeError, NotImplementedError) as exc:
            raise ExactComputationUnavailable("CAD could not construct a squarefree basis") from exc
        if not factors and primitive.total_degree() > 0:
            raise ExactComputationUnavailable(
                "CAD factorization lost a nonconstant polynomial"
            ) from None
        pieces = []
        for factor, _multiplicity in factors:
            try:
                piece = sp.Poly(factor, *primitive.gens, domain=primitive.domain)
            except (sp.PolynomialError, TypeError, ValueError):
                piece = sp.Poly(factor, *primitive.gens, extension=True)
            if piece.total_degree() > 0:
                pieces.append(piece)
        if not pieces:
            return None
        result = pieces[0]
        for piece in pieces[1:]:
            result *= piece

    result = result.primitive()[1]
    if result.LC().could_extract_minus_sign():
        result = -result
    return result


def lower_polynomial(expr: sp.Expr, lower_gens: Sequence[sp.Symbol]) -> sp.Poly | None:
    """Convert a projected expression to a normalized lower-dimensional polynomial."""

    if not lower_gens:
        return None
    try:
        poly = sp.Poly(expr, *lower_gens)
    except (sp.PolynomialError, ValueError, TypeError):
        return None
    return normalize_poly(poly)


def univariate_coefficients(poly: sp.Poly, var: sp.Symbol) -> tuple[sp.Expr, ...]:
    """Return coefficients when ``poly`` is viewed as univariate in ``var``."""

    return tuple(sp.Poly(poly.as_expr(), var).all_coeffs())


def univariate_leading_coefficient(poly: sp.Poly, var: sp.Symbol) -> sp.Expr:
    """Return the leading coefficient when viewed as univariate in ``var``."""

    return sp.Poly(poly.as_expr(), var).LC()


@lru_cache(maxsize=512)
def _discriminant_expression_cached(expr: sp.Expr, var: sp.Symbol) -> sp.Expr:
    return sp.discriminant(expr, var)


def discriminant_expression(poly: sp.Poly, var: sp.Symbol) -> sp.Expr:
    """Return the discriminant expression with respect to ``var``."""

    return _discriminant_expression_cached(poly.as_expr(), var)


@lru_cache(maxsize=512)
def _resultant_expression_cached(
    left_expr: sp.Expr, right_expr: sp.Expr, var: sp.Symbol
) -> sp.Expr:
    from ..algebraic.modular import modular_resultant_qq

    accelerated = modular_resultant_qq(left_expr, right_expr, var)
    if accelerated is not None:
        return accelerated.resultant
    return sp.resultant(left_expr, right_expr, var)


def resultant_expression(left: sp.Poly, right: sp.Poly, var: sp.Symbol) -> sp.Expr:
    """Return the exact resultant, preferring centralized FLINT arithmetic."""

    accelerated = flint_resultant(left, right, var)
    if accelerated is not None:
        return accelerated
    return _resultant_expression_cached(left.as_expr(), right.as_expr(), var)


def derivative_resultant_expression(poly: sp.Poly, var: sp.Symbol) -> sp.Expr:
    """Return ``Res(poly, d(poly)/dvar, var)`` without repeated caller conversions."""

    differentiate_generator = var in poly.gens
    if differentiate_generator:
        derivative = poly.diff(var)
    else:
        expr = poly.as_expr()
        derivative_expr = sp.diff(expr, var)
        derivative = sp.Poly(derivative_expr, *poly.gens, domain=poly.domain)
    accelerated = flint_resultant(poly, derivative, var)
    if accelerated is not None:
        return accelerated
    if differentiate_generator:
        expr = poly.as_expr()
        derivative_expr = derivative.as_expr()
    return sp.resultant(expr, derivative_expr, var)
