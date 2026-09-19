from dataclasses import replace

import sympy as sp

from semialg.algebraic.gtz import (
    certified_independent_localization,
    contract_localized_ideal,
    saturation_stabilization,
    verify_independent_localization_certificate,
    verify_localization_contraction_certificate,
    verify_saturation_stabilization_certificate,
)


def test_independent_localization_certifies_zero_dimensional_extension():
    x, y, z = sp.symbols("x y z")
    result = certified_independent_localization((x * y - z, y**2 - x), (x, y, z))
    assert len(result.independent_variables) == 1
    assert set(result.independent_variables + result.dependent_variables) == {x, y, z}
    assert result.quotient_dimension > 0
    assert verify_independent_localization_certificate(result.certificate)


def test_independent_localization_handles_hypersurface_component():
    x, y = sp.symbols("x y")
    result = certified_independent_localization((x * y,), (x, y))
    assert len(result.independent_variables) == 1
    assert len(result.dependent_variables) == 1
    assert result.localized_basis == (result.dependent_variables[0],)
    assert result.quotient_dimension == 1


def test_independent_localization_certificate_rejects_tampering():
    x, y = sp.symbols("x y")
    result = certified_independent_localization((x * y,), (x, y))
    bad = replace(result.certificate, quotient_dimension=2)
    assert not verify_independent_localization_certificate(bad)


def test_contraction_is_invariant_under_localized_presentation():
    u, x = sp.symbols("u x")
    result = contract_localized_ideal((u * x,), (u,), (x,))
    assert result.generators == (x,)
    assert result.certificate.canonical_localized_basis == (x,)
    assert verify_localization_contraction_certificate(result.certificate)


def test_contraction_clears_denominators_and_saturates():
    u, x = sp.symbols("u x")
    result = contract_localized_ideal((x + 1 / u,), (u,), (x,))
    assert result.generators == (u * x + 1,)
    assert sp.expand(result.denominator_product - u) == 0
    assert verify_localization_contraction_certificate(result.certificate)


def test_contraction_certificate_rejects_bad_denominator():
    u, x = sp.symbols("u x")
    result = contract_localized_ideal((x + 1 / u,), (u,), (x,))
    bad = replace(result.certificate, denominator_product=u + 1)
    assert not verify_localization_contraction_certificate(bad)


def test_saturation_stabilization_finds_exact_exponent_two():
    x, y = sp.symbols("x y")
    result = saturation_stabilization((x**2, x * y), x, (x, y))
    assert result.exponent == 2
    assert result.generators == (1,)
    assert result.companion_generators == (x**2, x * y)
    assert verify_saturation_stabilization_certificate(result.certificate)


def test_saturation_stabilization_recovers_gtz_split_identity():
    x, y = sp.symbols("x y")
    result = saturation_stabilization((x * y,), x, (x, y))
    assert result.exponent == 1
    assert result.generators == (y,)
    assert result.companion_generators == (x,)
    assert verify_saturation_stabilization_certificate(result.certificate)


def test_already_saturated_ideal_has_exponent_zero():
    x, y = sp.symbols("x y")
    result = saturation_stabilization((x * y,), x + y, (x, y))
    assert result.exponent == 0
    assert result.generators == (x * y,)
    assert result.companion_generators == (1,)


def test_saturation_certificate_rejects_tampered_exponent():
    x, y = sp.symbols("x y")
    result = saturation_stabilization((x * y,), x, (x, y))
    bad = replace(result.certificate, exponent=0)
    assert not verify_saturation_stabilization_certificate(bad)


def test_unified_replay_supports_gtz_results():
    from semialg import replay_certificate

    x, y = sp.symbols("x y")
    localization = certified_independent_localization((x * y,), (x, y))
    contraction = contract_localized_ideal((x,), (y,), (x,))
    stabilization = saturation_stabilization((x * y,), x, (x, y))
    assert replay_certificate(localization).verified is True
    assert replay_certificate(contraction).verified is True
    assert replay_certificate(stabilization).verified is True
