"""Systematic mutation tests for replayable certificate dataclasses."""

import sympy as sp

from semialg.algebraic.gtz import (
    certified_independent_localization,
    contract_localized_ideal,
    saturation_stabilization,
    verify_independent_localization_certificate,
    verify_localization_contraction_certificate,
    verify_saturation_stabilization_certificate,
)
from semialg.algebraic.modular import (
    modular_groebner_basis_qq,
    modular_resultant_qq,
    modular_subresultants_qq,
    verify_modular_groebner_certificate,
    verify_modular_resultant_certificate,
    verify_modular_subresultant_certificate,
)
from semialg.sos_certificates import SOSCertificate, verify_sos_certificate

from ._certificate_mutation import assert_certificate_fields_reject_mutation


def test_modular_groebner_certificate_proof_fields_reject_mutation():
    x, y = sp.symbols("x y")
    result = modular_groebner_basis_qq((x**2 - y, x * y - 1), (x, y))
    assert result is not None
    assert_certificate_fields_reject_mutation(
        result.certificate,
        verify_modular_groebner_certificate,
        ignored=frozenset({"certified"}),
    )


def test_modular_resultant_certificate_proof_fields_reject_mutation():
    a, x = sp.symbols("a x")
    result = modular_resultant_qq(x**2 + a, x + 1, x)
    assert result is not None
    assert_certificate_fields_reject_mutation(
        result.certificate,
        verify_modular_resultant_certificate,
        ignored=frozenset({"certified"}),
    )


def test_modular_subresultant_certificate_proof_fields_reject_mutation():
    a, x = sp.symbols("a x")
    result = modular_subresultants_qq(x**2 + a * x + 1, x + a, x)
    assert result is not None
    assert_certificate_fields_reject_mutation(
        result.certificate,
        verify_modular_subresultant_certificate,
        ignored=frozenset({"certified"}),
    )


def test_gtz_localization_certificate_proof_fields_reject_mutation():
    x, y = sp.symbols("x y")
    result = certified_independent_localization((x * y,), (x, y))
    assert_certificate_fields_reject_mutation(
        result.certificate, verify_independent_localization_certificate
    )


def test_gtz_contraction_certificate_proof_fields_reject_mutation():
    u, x = sp.symbols("u x")
    result = contract_localized_ideal((x + 1 / u,), (u,), (x,))
    assert_certificate_fields_reject_mutation(
        result.certificate, verify_localization_contraction_certificate
    )


def test_gtz_saturation_certificate_proof_fields_reject_mutation():
    x, y = sp.symbols("x y")
    result = saturation_stabilization((x * y,), x, (x, y))
    assert_certificate_fields_reject_mutation(
        result.certificate, verify_saturation_stabilization_certificate
    )


def test_sos_certificate_proof_fields_reject_mutation():
    x = sp.symbols("x")
    cert = SOSCertificate(x**2 + 1, (x,), (sp.Integer(1), x), sp.ImmutableMatrix.eye(2))
    assert_certificate_fields_reject_mutation(
        cert, lambda candidate: verify_sos_certificate(x**2 + 1, candidate)
    )


def test_regular_chain_certificate_proof_fields_reject_mutation():
    from semialg.algebraic_decomposition import (
        recursive_regular_chain_decomposition,
        verify_regular_chain_decomposition_certificate,
    )

    x, y = sp.symbols("x y")
    result = recursive_regular_chain_decomposition((x * y,), (x, y))
    assert result.certificate is not None
    assert_certificate_fields_reject_mutation(
        result.certificate,
        verify_regular_chain_decomposition_certificate,
        ignored=frozenset(),
    )


def test_radical_certificate_proof_fields_reject_mutation():
    from semialg.algebraic_decomposition import radical_ideal, verify_radical_ideal_certificate

    x, y = sp.symbols("x y")
    result = radical_ideal((x**2, x * y), (x, y))
    assert result.certificate is not None
    assert_certificate_fields_reject_mutation(result.certificate, verify_radical_ideal_certificate)


def test_minimal_prime_certificate_proof_fields_reject_mutation():
    from semialg.algebraic_decomposition import (
        certified_radical_minimal_prime_decomposition,
        verify_minimal_prime_decomposition_certificate,
    )

    x, y = sp.symbols("x y")
    result = certified_radical_minimal_prime_decomposition((x * y,), (x, y))
    assert result.certificate is not None
    assert_certificate_fields_reject_mutation(
        result.certificate,
        verify_minimal_prime_decomposition_certificate,
        ignored=frozenset(),
    )


def test_primary_certificate_proof_fields_reject_mutation():
    from semialg.algebraic_decomposition import (
        primary_decomposition,
        verify_primary_decomposition_certificate,
    )

    x = sp.symbols("x")
    result = primary_decomposition((x**2 * (x - 1),), (x,))
    assert result.certificate is not None
    assert_certificate_fields_reject_mutation(
        result.certificate,
        verify_primary_decomposition_certificate,
        ignored=frozenset(),
    )


def test_zero_dimensional_primary_certificate_proof_fields_reject_mutation():
    from semialg.algebraic.gtz_zero_dim import (
        verify_zero_dimensional_primary_certificate,
        zero_dimensional_primary_decomposition,
    )

    x = sp.symbols("x")
    result = zero_dimensional_primary_decomposition(((x - 1) ** 2 * (x + 2),), (x,))
    assert result.certificate is not None
    assert_certificate_fields_reject_mutation(
        result.certificate,
        verify_zero_dimensional_primary_certificate,
        ignored=frozenset({"coefficient_field", "field_relations"}),
    )


def test_recursive_gtz_certificate_proof_fields_reject_mutation():
    from semialg.algebraic.gtz_primary import (
        gtz_primary_decomposition,
        verify_gtz_primary_decomposition_certificate,
    )

    x, y = sp.symbols("x y")
    result = gtz_primary_decomposition((x * y,), (x, y))
    assert result.certificate is not None
    assert_certificate_fields_reject_mutation(
        result.certificate,
        verify_gtz_primary_decomposition_certificate,
        ignored=frozenset(),
    )
