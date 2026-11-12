"""Finite primary decomposition over exact coefficient fields and their towers."""

from dataclasses import dataclass

import sympy as sp
from algroots.basis_certificates import certify_groebner_basis
from algroots.domains import coefficient_domain, polynomial
from algroots.finite_decomposition import primary_decomposition as finite_primary
from algroots.quotient import QuotientAlgebra

from .._algebraic_decomposition_types import (
    PrimaryComponent,
    PrimaryDecompositionCertificate,
    PrimaryDecompositionResult,
)


@dataclass(frozen=True)
class CoefficientPrimaryCertificate(PrimaryDecompositionCertificate):
    """Replay payload for finite primary factors over a declared exact field."""

    coefficient_field: object
    basis_certificate: object

    def __reduce__(self):
        return _restore_certificate, (
            self.variables,
            self.source_equations,
            self.components,
            self.complete,
            self.irredundant,
            self.basis_certificate,
        )


def _restore_certificate(variables, source, components, complete, irredundant, proof):
    return CoefficientPrimaryCertificate(
        variables, source, components, complete, irredundant, proof.basis.domain, proof
    )


def coefficient_primary_decomposition(
    generators, variables, field, *, backend="native", limits=None
):
    """Compute finite primary factors in the declared field without a QQ lift.

    Rational parameter fields and monogenic towers are normalized through their
    exact domain interface. An explicit native choice applies to the source
    basis; component construction and replay use exact coefficient arithmetic.
    """
    from algroots._validation import normalize_equation, validate_variables

    domain = coefficient_domain(field)
    variables = validate_variables(variables)
    source = tuple(
        polynomial(normalize_equation(p), variables, domain).as_expr() for p in generators
    )
    algebra = QuotientAlgebra.from_polynomials(
        source, variables, domain=domain, groebner_backend=backend, groebner_limits=limits
    )
    domain = algebra.domain
    proof = certify_groebner_basis(
        source, variables, candidate=algebra.groebner_basis, domain=domain
    )
    components = tuple(
        PrimaryComponent(
            tuple(p.as_expr() for p in component.groebner_basis.polys),
            tuple(p.as_expr() for p in residue.groebner_basis.polys),
            0,
            component.dimension,
            True,
            "coefficient_field_primary",
        )
        for component, residue in finite_primary(
            algebra, groebner_backend=backend, groebner_limits=limits
        )
    )
    certificate = CoefficientPrimaryCertificate(
        variables, source, components, True, True, domain, proof
    )
    result = PrimaryDecompositionResult(
        variables, components, True, True, "coefficient_field_primary", certificate
    )
    if not verify_coefficient_primary_certificate(certificate):
        raise ArithmeticError("coefficient-field primary decomposition failed exact replay")
    return result


def verify_coefficient_primary_certificate(certificate):
    """Recheck source ideal, primary factor theorem and every returned component."""
    try:
        if (
            not isinstance(certificate, CoefficientPrimaryCertificate)
            or not certificate.complete
            or not certificate.irredundant
        ):
            return False
        proof = certificate.basis_certificate
        if (
            not proof.verify()
            or proof.basis.domain != certificate.coefficient_field
            or proof.basis.gens != certificate.variables
        ):
            return False
        domain = certificate.coefficient_field
        source = tuple(
            polynomial(p, certificate.variables, domain) for p in certificate.source_equations
        )
        if source != tuple(polynomial(p, certificate.variables, domain) for p in proof.source):
            return False
        algebra = QuotientAlgebra.from_certificate(proof)
        expected = finite_primary(algebra, groebner_backend="sympy")
        if len(expected) != len(certificate.components):
            return False
        for component, (factor, residue) in zip(certificate.components, expected, strict=True):
            if (
                component.dimension != 0
                or component.degree != factor.dimension
                or not component.primary
                or component.method != "coefficient_field_primary"
            ):
                return False
            actual = QuotientAlgebra.from_polynomials(
                component.equations, certificate.variables, domain=domain, groebner_backend="sympy"
            )
            radical = QuotientAlgebra.from_polynomials(
                component.radical, certificate.variables, domain=domain, groebner_backend="sympy"
            )
            if (
                actual.groebner_basis != factor.groebner_basis
                or radical.groebner_basis != residue.groebner_basis
            ):
                return False
        return True
    except (
        ValueError,
        TypeError,
        AttributeError,
        ArithmeticError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False
