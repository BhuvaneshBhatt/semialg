from dataclasses import replace

import sympy as sp

from semialg.algebraic_decomposition import (
    certified_radical_minimal_prime_decomposition,
    radical_ideal,
    recursive_regular_chain_decomposition,
    verify_minimal_prime_decomposition_certificate,
    verify_radical_ideal_certificate,
    verify_regular_chain_decomposition_certificate,
)


def test_regular_chain_certificate_replays_without_search():
    x, y, z = sp.symbols("x y z")
    result = recursive_regular_chain_decomposition((x * y + z, y * z), (x, y, z))
    assert result.certificate is not None
    assert result.certificate.splits
    assert verify_regular_chain_decomposition_certificate(result.certificate)


def test_regular_chain_certificate_rejects_tampered_child():
    x, y = sp.symbols("x y")
    result = recursive_regular_chain_decomposition((x * y,), (x, y))
    certificate = result.certificate
    assert certificate is not None and certificate.splits
    branch = certificate.splits[0]
    bad_branch = replace(branch, children=(branch.children[0], (x + y,)))
    bad = replace(certificate, splits=(bad_branch, *certificate.splits[1:]))
    assert not verify_regular_chain_decomposition_certificate(bad)


def test_regular_chain_certificate_rejects_tampered_degree():
    x, y = sp.symbols("x y")
    result = recursive_regular_chain_decomposition((x * y,), (x, y))
    certificate = result.certificate
    assert certificate is not None
    bad_component = replace(certificate.components[0], degree=99)
    bad = replace(certificate, components=(bad_component, *certificate.components[1:]))
    assert not verify_regular_chain_decomposition_certificate(bad)


def test_radical_certificate_replays_and_rejects_generator_tamper():
    x, y = sp.symbols("x y")
    result = radical_ideal((x**2, x * y), (x, y))
    assert result.certificate is not None
    assert verify_radical_ideal_certificate(result.certificate)
    bad = replace(result.certificate, generators=(x * y,))
    assert not verify_radical_ideal_certificate(bad)


def test_minimal_prime_certificate_replays_and_rejects_prime_promotion():
    x, y = sp.symbols("x y")
    result = certified_radical_minimal_prime_decomposition((x * y,), (x, y))
    assert result.certificate is not None
    assert verify_minimal_prime_decomposition_certificate(result.certificate)
    bad_component = replace(result.certificate.components[0], degree=7)
    bad = replace(
        result.certificate,
        components=(bad_component, *result.certificate.components[1:]),
    )
    assert not verify_minimal_prime_decomposition_certificate(bad)


def test_unified_replay_certificate_uses_algebraic_verifier():
    from semialg.certificates import replay_certificate

    x, y = sp.symbols("x y")
    result = radical_ideal((x**2, x * y), (x, y))
    replay = replay_certificate(result)
    assert replay.verified is True
    assert replay.kind == "radical_ideal"


def test_regular_chain_certificate_rejects_tampered_saturation_splitter():
    x, y, z = sp.symbols("x y z")
    result = recursive_regular_chain_decomposition((x * y + z, y * z), (x, y, z))
    certificate = result.certificate
    assert certificate is not None and certificate.splits[0].splitter is not None
    bad_branch = replace(certificate.splits[0], splitter=x + y + 1)
    bad = replace(certificate, splits=(bad_branch, *certificate.splits[1:]))
    assert not verify_regular_chain_decomposition_certificate(bad)
