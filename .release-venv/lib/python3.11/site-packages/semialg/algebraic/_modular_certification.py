"""Independent replay of certificates produced by modular algebra algorithms."""

from __future__ import annotations

import sympy as sp

from ._modular_resultants import (
    _input_expr,
    _parameter_degree_bounds,
    _resultant_grid,
    _sylvester_det_at,
)
from ._modular_types import (
    ModularFractionFieldGroebnerCertificate,
    ModularGroebnerCertificate,
    ModularResultantCertificate,
    ModularSubresultantCertificate,
)
from .modular import _canonical_exact_basis


def verify_modular_groebner_certificate(
    certificate: ModularGroebnerCertificate,
) -> bool:
    """Replay the exact proof obligations of a modular Groebner certificate."""
    if not certificate.certified:
        return False
    try:
        gb, canonical = _canonical_exact_basis(
            certificate.basis, certificate.variables, certificate.order
        )
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        NotImplementedError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False
    normalized = tuple(sp.expand(p) for p in certificate.basis)
    if set(canonical) != set(normalized) or len(canonical) != len(normalized):
        return False
    for generator in certificate.source_generators:
        try:
            _, remainder = gb.reduce(generator)
        except (
            ArithmeticError,
            ValueError,
            TypeError,
            NotImplementedError,
            sp.PolynomialError,
            sp.polys.polyerrors.CoercionFailed,
        ):
            return False
        if sp.expand(remainder) != 0:
            return False
    if len(certificate.membership_representations) != len(certificate.basis):
        return False
    for polynomial, multipliers in zip(
        certificate.basis, certificate.membership_representations, strict=False
    ):
        if len(multipliers) != len(certificate.source_generators):
            return False
        rebuilt = sum(
            h * g for h, g in zip(multipliers, certificate.source_generators, strict=False)
        )
        if sp.expand(polynomial - rebuilt) != 0:
            return False
    return True


def verify_modular_fraction_field_groebner_certificate(
    certificate: ModularFractionFieldGroebnerCertificate,
) -> bool:
    """Replay only exact QQ(U) proof obligations; modular images are advisory."""
    if not certificate.certified:
        return False
    domain = sp.QQ.frac_field(*certificate.parameters)
    try:
        gb = sp.groebner(
            certificate.basis,
            *certificate.variables,
            order=certificate.order,
            domain=domain,
        )
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False
    canonical = tuple(sp.cancel(poly.as_expr()) for poly in gb.polys)
    if canonical != tuple(sp.cancel(poly) for poly in certificate.basis):
        return False
    for generator in certificate.source_generators:
        try:
            _quotients, remainder = gb.reduce(generator)
        except (
            ArithmeticError,
            ValueError,
            TypeError,
            sp.PolynomialError,
            sp.polys.polyerrors.CoercionFailed,
        ):
            return False
        if sp.cancel(remainder) != 0:
            return False
    if len(certificate.membership_representations) != len(certificate.basis):
        return False
    for polynomial, multipliers in zip(
        certificate.basis, certificate.membership_representations, strict=True
    ):
        if len(multipliers) != len(certificate.source_generators):
            return False
        rebuilt = sum(
            multiplier * generator
            for multiplier, generator in zip(
                multipliers, certificate.source_generators, strict=True
            )
        )
        if sp.cancel(polynomial - rebuilt) != 0:
            return False
    return True


def verify_modular_resultant_certificate(certificate: ModularResultantCertificate) -> bool:
    """Replay a resultant certificate without symbolic resultant computation."""
    if not certificate.certified:
        return False
    try:
        all_vars = (certificate.variable, *certificate.parameters)
        f = sp.Poly(certificate.first, *all_vars, domain=sp.QQ)
        g = sp.Poly(certificate.second, *all_vars, domain=sp.QQ)
        candidate = (
            sp.Poly(certificate.resultant, *certificate.parameters, domain=sp.QQ)
            if certificate.parameters
            else None
        )
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        NotImplementedError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False
    m = f.degree(certificate.variable)
    n = g.degree(certificate.variable)
    expected_bounds = _parameter_degree_bounds(
        certificate.first, certificate.second, certificate.variable, certificate.parameters
    )
    if tuple(certificate.parameter_degree_bounds) != expected_bounds:
        return False
    if certificate.parameters:
        for parameter, bound in zip(certificate.parameters, expected_bounds, strict=False):
            if candidate.degree(parameter) > bound:
                return False
    expected_grid = _resultant_grid(certificate.parameters, expected_bounds)
    if tuple(certificate.verification_grid) != expected_grid:
        return False
    for point in expected_grid:
        substitutions = dict(zip(certificate.parameters, point, strict=False))
        expected = _sylvester_det_at(
            certificate.first, certificate.second, certificate.variable, m, n, substitutions
        )
        actual = sp.factor(_input_expr(certificate.resultant).subs(substitutions))
        if sp.expand(actual - expected) != 0:
            return False
    return True


def verify_modular_subresultant_certificate(
    certificate: ModularSubresultantCertificate,
) -> bool:
    """Replay exact normalization and resultant obligations for a PRS certificate."""
    if not certificate.certified:
        return False
    if not verify_modular_resultant_certificate(certificate.resultant_certificate):
        return False
    try:
        domain = sp.QQ.poly_ring(*certificate.parameters) if certificate.parameters else sp.QQ
        left = sp.Poly(certificate.first, certificate.variable, domain=domain)
        right = sp.Poly(certificate.second, certificate.variable, domain=domain)
        if left.degree() < right.degree():
            left, right = right, left
        exact = tuple(sp.subresultants(left, right))
        exact_exprs = tuple(sp.expand(p.as_expr() if isinstance(p, sp.Poly) else p) for p in exact)
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        NotImplementedError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False
    candidates = tuple(sp.expand(p) for p in certificate.polynomials)
    return candidates == exact_exprs
