"""Certified zero-dimensional primary decomposition over algebraic fields.

For a finite algebraic extension ``K/QQ`` this module lifts an ideal in
``K[X]`` to ``QQ[alpha_1,...,alpha_r,X]`` by adjoining the defining equations
of the coefficient-field tower.  Zero-dimensional primary components are then
isolated by certified maximal-prime separators and exact saturation.  The
components are transported back to ``K[X]`` and verified by lifting again.

The primary theorem used here is specific and strong: in a zero-dimensional
Noetherian quotient all prime ideals are maximal, hence there are no embedded
prime inclusions; saturating away all other maximal primes isolates a primary
component for the remaining maximal prime.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ..algebraic_function_fields import (
    MonogenicFunctionField,
    RationalFunctionField,
    certified_factor_univariate,
    field_element_to_expr,
    tower_degree,
    tower_evaluate,
)
from .equality_ideal import EqualityIdealContext
from .gtz import (
    SaturationStabilizationCertificate,
    saturation_stabilization,
    verify_saturation_stabilization_certificate,
)
from .hilbert import ideal_degree
from .ideal_ops import canonical_qq_basis_uncached, intersection_all_qq, qq_ideal_equal_uncached


def _canonical_qq_basis(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    return canonical_qq_basis_uncached(generators, variables)


def _qq_ideal_equal(
    left: Sequence[sp.Expr], right: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> bool:
    return qq_ideal_equal_uncached(left, right, variables)


def _intersection_all(
    ideals: Sequence[Sequence[sp.Expr]], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    return intersection_all_qq(ideals, variables)


def _tower_levels(field) -> tuple[MonogenicFunctionField, ...]:
    levels: list[MonogenicFunctionField] = []
    current = field
    while isinstance(current, MonogenicFunctionField):
        levels.append(current)
        current = current.base
    if not isinstance(current, RationalFunctionField) or current.parameters:
        raise ValueError(
            "zero-dimensional decomposition supports QQ and finite algebraic extensions of QQ"
        )
    levels.reverse()
    return tuple(levels)


def _coefficient_generators(field) -> tuple[sp.Symbol, ...]:
    return tuple(level.generator for level in _tower_levels(field))


def _field_relations(field) -> tuple[sp.Expr, ...]:
    relations = []
    for level in _tower_levels(field):
        relation = sum(
            field_element_to_expr(coefficient, level.base) * level.generator**exponent
            for exponent, coefficient in enumerate(level.modulus)
        )
        relations.append(sp.expand(relation))
    return tuple(relations)


def _verify_coefficient_field(field) -> bool:
    """Certify every defining polynomial irreducible over the preceding field."""
    try:
        for level in _tower_levels(field):
            if not certified_factor_univariate(level.base, level.modulus).irreducible:
                return False
        return True
    except (ArithmeticError, TypeError, ValueError, ZeroDivisionError):
        return False


def _canonical_field_polynomial(
    expression: sp.Expr, variables: Sequence[sp.Symbol], field
) -> sp.Expr:
    vars_ = tuple(variables)
    expression = sp.expand(sp.sympify(expression))
    allowed = set(vars_) | set(_coefficient_generators(field))
    if not expression.free_symbols.issubset(allowed):
        raise ValueError("polynomial contains symbols outside variables/coefficient field")
    if not vars_:
        return sp.expand(field_element_to_expr(tower_evaluate(expression, field), field))
    poly = sp.Poly(expression, *vars_, domain=sp.EX)
    out = sp.Integer(0)
    for monomial, coefficient in poly.terms():
        coefficient_expr = field_element_to_expr(tower_evaluate(coefficient, field), field)
        term = coefficient_expr
        for variable, exponent in zip(vars_, monomial, strict=True):
            term *= variable ** int(exponent)
        out += term
    return sp.expand(out)


def _canonical_field_generators(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol], field
) -> tuple[sp.Expr, ...]:
    out: list[sp.Expr] = []
    seen: set[sp.Expr] = set()
    for generator in generators:
        reduced = _canonical_field_polynomial(generator, variables, field)
        if reduced == 0 or reduced in seen:
            continue
        seen.add(reduced)
        out.append(reduced)
    return tuple(out)


def _separator_for_prime(index: int, primes: Sequence[EqualityIdealContext]) -> sp.Expr | None:
    target = primes[index]
    factors = []
    for j, other in enumerate(primes):
        if j == index:
            continue
        chosen = next((g for g in other.generators if target.normal_form(g) != 0), None)
        if chosen is None:
            return None
        factors.append(chosen)
    return sp.expand(sp.prod(factors)) if factors else sp.Integer(1)


@dataclass(frozen=True)
class ZeroDimensionalPrimaryComponent:
    """One certified primary component over the requested coefficient field."""

    generators: tuple[sp.Expr, ...]
    radical: tuple[sp.Expr, ...]
    degree: int
    primary: bool
    method: str
    ambient_generators: tuple[sp.Expr, ...]
    ambient_radical: tuple[sp.Expr, ...]
    separator: sp.Expr
    saturation_certificate: SaturationStabilizationCertificate


@dataclass(frozen=True)
class ZeroDimensionalPrimaryCertificate:
    """Replay payload for zero-dimensional primary decomposition."""

    variables: tuple[sp.Symbol, ...]
    coefficient_field: object
    source_generators: tuple[sp.Expr, ...]
    field_relations: tuple[sp.Expr, ...]
    ambient_variables: tuple[sp.Symbol, ...]
    components: tuple[ZeroDimensionalPrimaryComponent, ...]
    minimal_prime_cert: object
    extension_degree: int


@dataclass(frozen=True)
class ZeroDimensionalPrimaryResult:
    """Certified primary decomposition over QQ or a finite algebraic extension."""

    variables: tuple[sp.Symbol, ...]
    coefficient_field: object
    components: tuple[ZeroDimensionalPrimaryComponent, ...]
    complete: bool
    irredundant: bool
    method: str
    certificate: ZeroDimensionalPrimaryCertificate | None = None


def zero_dimensional_primary_decomposition(
    generators: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    coefficient_field: object | None = None,
) -> ZeroDimensionalPrimaryResult:
    """Compute primary components for a zero-dimensional ideal in ``K[X]``."""
    from ..algebraic_decomposition import certified_radical_minimal_prime_decomposition

    vars_ = tuple(variables)
    field = coefficient_field or RationalFunctionField(())
    coefficient_symbols = _coefficient_generators(field)
    if set(vars_).intersection(coefficient_symbols):
        raise ValueError("polynomial variables must be disjoint from coefficient generators")
    if not _verify_coefficient_field(field):
        raise ValueError("coefficient tower was not certified as a finite algebraic field")

    source = _canonical_field_generators(tuple(generators), vars_, field)
    relations = _field_relations(field)
    ambient_variables = (*coefficient_symbols, *vars_)
    ambient_source = (*relations, *source)
    context = EqualityIdealContext.build(ambient_source, ambient_variables)
    if context.inconsistent:
        # Unit ideal: empty primary decomposition.
        return ZeroDimensionalPrimaryResult(
            vars_,
            field,
            tuple(),
            True,
            True,
            "gtz3_algebraic_lift" if coefficient_symbols else "gtz3_qq",
            None,
        )
    if context.dimension != 0:
        raise ValueError(
            "zero-dimensional primary decomposition requires a zero-dimensional ideal over the coefficient field"
        )

    minimal = certified_radical_minimal_prime_decomposition(ambient_source, ambient_variables)
    if not minimal.minimal_primes_complete or minimal.certificate is None:
        return ZeroDimensionalPrimaryResult(
            vars_, field, tuple(), False, False, "minimal_prime_decomposition_incomplete"
        )
    primes = [
        EqualityIdealContext.build(component.equations, ambient_variables)
        for component in minimal.components
    ]
    if any(prime.dimension != 0 for prime in primes):
        raise ArithmeticError("zero-dimensional source produced a positive-dimensional prime")

    extension_degree = tower_degree(field)
    components: list[ZeroDimensionalPrimaryComponent] = []
    for index, prime in enumerate(primes):
        separator = _separator_for_prime(index, primes)
        if separator is None:
            return ZeroDimensionalPrimaryResult(
                vars_, field, tuple(), False, False, "separator_construction_incomplete"
            )
        saturation = saturation_stabilization(ambient_source, separator, ambient_variables)
        ambient_generators = tuple(saturation.generators)
        ambient_radical = tuple(prime.generators)
        ambient_degree = ideal_degree(ambient_generators, ambient_variables)
        if ambient_degree % extension_degree:
            raise ArithmeticError("component degree is incompatible with coefficient-field degree")
        mapped = _canonical_field_generators(ambient_generators, vars_, field)
        mapped_radical = _canonical_field_generators(ambient_radical, vars_, field)
        components.append(
            ZeroDimensionalPrimaryComponent(
                mapped,
                mapped_radical,
                ambient_degree // extension_degree,
                True,
                "zero_dimensional_maximal_localization",
                ambient_generators,
                ambient_radical,
                separator,
                saturation.certificate,
            )
        )

    if not _qq_ideal_equal(
        _intersection_all(
            [component.ambient_generators for component in components],
            ambient_variables,
        ),
        ambient_source,
        ambient_variables,
    ):
        raise ArithmeticError("primary components failed exact source reconstruction")

    certificate = ZeroDimensionalPrimaryCertificate(
        vars_,
        field,
        source,
        relations,
        ambient_variables,
        tuple(components),
        minimal.certificate,
        extension_degree,
    )
    result = ZeroDimensionalPrimaryResult(
        vars_,
        field,
        tuple(components),
        True,
        True,
        "gtz3_algebraic_lift" if coefficient_symbols else "gtz3_qq",
        certificate,
    )
    if not verify_zero_dimensional_primary_certificate(certificate):
        raise ArithmeticError("internal zero-dimensional primary certificate verification failed")
    return result


def verify_zero_dimensional_primary_certificate(
    certificate: ZeroDimensionalPrimaryCertificate,
) -> bool:
    """Replay field validity, maximal-prime separation, saturation, and transport."""
    from ..algebraic_decomposition import verify_minimal_prime_decomposition_certificate

    try:
        vars_ = tuple(certificate.variables)
        field = certificate.coefficient_field
        if not _verify_coefficient_field(field):
            return False
        coefficient_symbols = _coefficient_generators(field)
        ambient_variables = (*coefficient_symbols, *vars_)
        if ambient_variables != tuple(certificate.ambient_variables):
            return False
        relations = _field_relations(field)
        if relations != tuple(certificate.field_relations):
            return False
        extension_degree = tower_degree(field)
        if extension_degree != certificate.extension_degree:
            return False
        source = _canonical_field_generators(certificate.source_generators, vars_, field)
        if source != tuple(certificate.source_generators):
            return False
        ambient_source = (*relations, *source)
        context = EqualityIdealContext.build(ambient_source, ambient_variables)
        if context.inconsistent or context.dimension != 0:
            return False

        minimal_certificate = certificate.minimal_prime_cert
        if not verify_minimal_prime_decomposition_certificate(minimal_certificate):
            return False
        if not minimal_certificate.minimal_primes_complete:
            return False
        if not _qq_ideal_equal(
            minimal_certificate.source_equations, ambient_source, ambient_variables
        ):
            return False
        if len(certificate.components) != len(minimal_certificate.components):
            return False
        primes = [
            EqualityIdealContext.build(component.equations, ambient_variables)
            for component in minimal_certificate.components
        ]
        if any(prime.dimension != 0 for prime in primes):
            return False

        for index, (component, prime_component) in enumerate(
            zip(certificate.components, minimal_certificate.components, strict=True)
        ):
            if not component.primary or component.method != "zero_dimensional_maximal_localization":
                return False
            if tuple(component.ambient_radical) != tuple(prime_component.equations):
                return False
            separator = component.separator
            target = primes[index]
            if target.normal_form(separator) == 0:
                return False
            for j, other in enumerate(primes):
                if j != index and other.normal_form(separator) != 0:
                    return False
            sat = component.saturation_certificate
            if tuple(sat.variables) != ambient_variables:
                return False
            if not _qq_ideal_equal(sat.source_generators, ambient_source, ambient_variables):
                return False
            if sp.expand(sat.splitter - separator) != 0:
                return False
            if not verify_saturation_stabilization_certificate(sat):
                return False
            if not _qq_ideal_equal(
                component.ambient_generators, sat.saturated_generators, ambient_variables
            ):
                return False
            mapped = _canonical_field_generators(component.ambient_generators, vars_, field)
            radical = _canonical_field_generators(component.ambient_radical, vars_, field)
            if mapped != tuple(component.generators) or radical != tuple(component.radical):
                return False
            if not _qq_ideal_equal(
                (*relations, *component.generators),
                component.ambient_generators,
                ambient_variables,
            ):
                return False
            if not _qq_ideal_equal(
                (*relations, *component.radical),
                component.ambient_radical,
                ambient_variables,
            ):
                return False
            ambient_degree = ideal_degree(component.ambient_generators, ambient_variables)
            if ambient_degree % extension_degree:
                return False
            if component.degree != ambient_degree // extension_degree:
                return False

        reconstructed = _intersection_all(
            [component.ambient_generators for component in certificate.components],
            ambient_variables,
        )
        return _qq_ideal_equal(reconstructed, ambient_source, ambient_variables)
    except (
        ArithmeticError,
        TypeError,
        ValueError,
        ZeroDivisionError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False


__all__ = [
    "ZeroDimensionalPrimaryCertificate",
    "ZeroDimensionalPrimaryComponent",
    "ZeroDimensionalPrimaryResult",
    "verify_zero_dimensional_primary_certificate",
    "zero_dimensional_primary_decomposition",
]
