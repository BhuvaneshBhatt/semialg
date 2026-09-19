"""Cross-cutting contracts for replayable exact certificates."""

from __future__ import annotations

import pickle
from dataclasses import replace

import pytest
import sympy as sp

import semialg
from semialg.algebraic import (
    certify_polynomial_root_interval,
    modular_groebner_basis_qq,
    verify_modular_groebner_certificate,
    verify_polynomial_root_interval_certificate,
)
from semialg.algebraic.gtz import (
    certified_independent_localization,
    saturation_stabilization,
    verify_independent_localization_certificate,
    verify_saturation_stabilization_certificate,
)
from semialg.cache_control import clear_caches
from semialg.polyhedral import HRepresentation, verify_h_redundancy_certificate

from ._certificate_mutation import assert_certificate_fields_reject_mutation


def _root_case():
    x = sp.Symbol("x")
    certificate = certify_polynomial_root_interval((x**2 - 2) * (x - 3), -2, 2, var=x)
    return certificate, verify_polynomial_root_interval_certificate


def _modular_case():
    x, y = sp.symbols("x y")
    result = modular_groebner_basis_qq((x**2 - y, x * y - 1), (x, y))
    assert result is not None
    return result.certificate, verify_modular_groebner_certificate


def _localization_case():
    x, y = sp.symbols("x y")
    result = certified_independent_localization((x * y,), (x, y))
    return result.certificate, verify_independent_localization_certificate


def _saturation_case():
    x, y = sp.symbols("x y")
    result = saturation_stabilization((x * y,), x, (x, y))
    return result.certificate, verify_saturation_stabilization_certificate


@pytest.mark.parametrize(
    "factory",
    (_root_case, _modular_case, _localization_case, _saturation_case),
    ids=("root-interval", "modular-groebner", "localization", "saturation"),
)
def test_replayable_certificates_survive_pickle_and_cache_reset(factory):
    certificate, verifier = factory()
    assert verifier(certificate)

    restored = pickle.loads(pickle.dumps(certificate))
    clear_caches(include_sympy=False, collect=False)

    assert type(restored) is type(certificate)
    assert verifier(restored)


def test_root_interval_certificate_rejects_every_proof_field_mutation():
    certificate, _ = _root_case()

    def verifier(candidate):
        return verify_polynomial_root_interval_certificate(
            candidate,
            polynomial=certificate.polynomial,
            var=certificate.polynomial.gens[0],
            left=certificate.left,
            right=certificate.right,
            include_left=certificate.include_left,
            include_right=certificate.include_right,
        )

    assert_certificate_fields_reject_mutation(certificate, verifier)


def test_root_interval_replay_binds_polynomial_variable_and_endpoint_semantics():
    x, y = sp.symbols("x y")
    certificate = certify_polynomial_root_interval(x**2 - 2, 0, 2, var=x)

    assert verify_polynomial_root_interval_certificate(certificate)
    assert not verify_polynomial_root_interval_certificate(certificate, var=y)
    assert not verify_polynomial_root_interval_certificate(certificate, polynomial=x**2 - 3, var=x)
    assert not verify_polynomial_root_interval_certificate(
        replace(certificate, include_left=not certificate.include_left),
        include_left=certificate.include_left,
    )


def test_h_redundancy_certificate_is_bound_to_its_representation():
    representation = HRepresentation(((-1,), (1,), (2,)), (0, 1, 2))
    certificates = representation.redundancy_certificates()
    certificate = certificates[2]
    assert verify_h_redundancy_certificate(representation, certificate)

    different = HRepresentation(((-1,), (1,), (2,)), (0, 1, 1))
    assert not verify_h_redundancy_certificate(different, certificate)


def test_every_public_certificate_type_has_an_explicit_replay_policy():
    import inspect

    import semialg.algebraic as algebraic
    import semialg.polyhedral as polyhedral
    import semialg.sos_certificates as sos_certificates

    exported = {
        name
        for module in (semialg, algebraic, polyhedral, sos_certificates)
        for name in module.__all__
        if name.endswith("Certificate") and inspect.isclass(getattr(module, name))
    }
    policies = {
        "HConstraintRedundancyCertificate": "standalone",
        "SOSCertificate": "standalone",
        "FractionFieldZeroDimensionalCertificate": "standalone",
        "GTZContractedComponentCertificate": "nested-gtz",
        "GTZNodeCertificate": "nested-gtz",
        "GTZPrimaryDecompositionCertificate": "standalone",
        "IndependentLocalizationCertificate": "standalone",
        "LocalizationContractionCertificate": "standalone",
        "ModularFractionFieldGroebnerCertificate": "standalone",
        "ModularGroebnerCertificate": "standalone",
        "ModularResultantCertificate": "standalone",
        "ModularSubresultantCertificate": "standalone",
        "PolynomialRootIntervalCertificate": "standalone",
        "SaturationStabilizationCertificate": "standalone",
        "ZeroDimensionalPrimaryCertificate": "standalone",
    }

    assert exported == policies.keys()
    assert set(policies.values()) <= {"standalone", "nested-gtz"}
