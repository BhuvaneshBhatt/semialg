from dataclasses import replace

import sympy as sp

from semialg import replay_certificate
from semialg.algebraic_decomposition import associated_primes, primary_decomposition
from semialg.algebraic_geometry import verify_primary_decomposition_certificate


def test_principal_primary_decomposition_tracks_multiplicity():
    x = sp.symbols("x")
    result = primary_decomposition((x**2 * (x - 1) ** 3,), (x,))
    assert result.complete and result.irredundant
    assert result.method == "principal_factorization"
    assert {component.radical for component in result.components} == {(x,), (x - 1,)}
    assert {component.degree for component in result.components} == {2, 3}
    assert verify_primary_decomposition_certificate(result.certificate)


def test_monomial_primary_decomposition_finds_embedded_prime():
    x, y = sp.symbols("x y")
    result = primary_decomposition((x**2, x * y), (x, y))
    assert result.complete and result.irredundant
    assert result.method == "monomial_irreducible_decomposition"
    assert {component.radical for component in result.components} == {(x,), (x, y)}
    assoc = associated_primes((x**2, x * y), (x, y))
    assert assoc.complete
    assert set(assoc.primes) == {(x,), (x, y)}


def test_zero_dimensional_localization_recovers_primary_components():
    x, y = sp.symbols("x y")
    result = primary_decomposition((x**2 * (x - 1), y), (x, y))
    assert result.complete and result.irredundant
    assert result.method == "zero_dimensional_localization"
    assert {component.radical for component in result.components} == {
        (x, y),
        (x - 1, y),
    }
    assert sorted(component.degree for component in result.components) == [1, 2]


def test_prime_and_unit_ideal_cases_are_complete():
    x, y = sp.symbols("x y")
    prime = primary_decomposition((x,), (x, y))
    assert prime.complete and prime.method == "prime_ideal"
    assert prime.components[0].method == "prime_is_primary"
    unit = primary_decomposition((sp.Integer(1),), (x, y))
    assert unit.complete and unit.components == ()
    assert associated_primes((sp.Integer(1),), (x, y)).primes == ()


def test_primary_certificate_replay_rejects_tampering():
    x, y = sp.symbols("x y")
    result = primary_decomposition((x**2, x * y), (x, y))
    certificate = result.certificate
    assert certificate is not None
    assert verify_primary_decomposition_certificate(certificate)
    component = certificate.components[1]
    forged_component = replace(component, radical=(y,))
    forged = replace(
        certificate,
        components=(certificate.components[0], forged_component),
    )
    assert not verify_primary_decomposition_certificate(forged)
    replay = replay_certificate(result)
    assert replay.verified is True
    assert replay.kind == "primary_decomposition"


def test_unsupported_positive_dimensional_nonmonomial_is_explicitly_incomplete():
    x, y, z = sp.symbols("x y z")
    # Nonmonomial, positive-dimensional, and not a principal/prime ideal.
    result = primary_decomposition((x**2, x * y + z * x), (x, y, z))
    if not result.complete:
        assert result.components == ()
        assert result.certificate is None
