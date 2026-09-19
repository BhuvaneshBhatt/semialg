import sympy as sp

from semialg.algebraic.groebner_utils import compute_groebner_basis
from semialg.algebraic.modular import (
    ModularGroebnerCertificate,
    modular_groebner_basis_qq,
    verify_modular_groebner_certificate,
)


def _exprs(basis):
    return tuple(sp.expand(poly.as_expr()) for poly in basis.polys)


def test_modular_groebner_reconstructs_and_certifies_grevlex():
    x, y = sp.symbols("x y")
    generators = (x * y - 1, y**2 - x)
    result = modular_groebner_basis_qq(generators, (x, y), order="grevlex")
    assert result is not None
    assert result.certified
    assert len(result.certificate.primes) >= 2
    assert verify_modular_groebner_certificate(result.certificate)

    direct = sp.groebner(generators, x, y, order="grevlex", domain=sp.QQ)
    assert result.basis == _exprs(direct)


def test_modular_groebner_handles_rational_coefficients_and_lex():
    x, y = sp.symbols("x y")
    generators = (sp.Rational(2, 3) * x**2 + y, x * y - sp.Rational(5, 7))
    result = modular_groebner_basis_qq(generators, (x, y), order="lex")
    assert result is not None
    assert verify_modular_groebner_certificate(result.certificate)
    direct = sp.groebner(generators, x, y, order="lex", domain=sp.QQ)
    assert result.basis == _exprs(direct)


def test_modular_groebner_certificate_rejects_tampered_membership():
    x, y = sp.symbols("x y")
    result = modular_groebner_basis_qq((x**2 - y, x * y - 1), (x, y))
    assert result is not None
    cert = result.certificate
    bad_reps = list(cert.membership_representations)
    bad_reps[0] = tuple(sp.Integer(0) for _ in cert.source_generators)
    tampered = ModularGroebnerCertificate(
        cert.variables,
        cert.order,
        cert.source_generators,
        cert.basis,
        tuple(bad_reps),
        cert.primes,
        cert.modulus,
    )
    assert not verify_modular_groebner_certificate(tampered)


def test_shared_groebner_helper_matches_direct_exact_basis():
    x, y = sp.symbols("x y")
    generators = (2 * x**2 + 3 * y, x * y - 5)
    accelerated = compute_groebner_basis(generators, (x, y), order="grevlex", domain=sp.QQ)
    direct = compute_groebner_basis(
        generators, (x, y), order="grevlex", domain=sp.QQ, modular=False
    )
    assert _exprs(accelerated) == _exprs(direct)


def test_shared_helper_falls_back_for_parameter_coefficients():
    a, x, y = sp.symbols("a x y")
    generators = (a * x + y, x**2 - 1)
    accelerated = compute_groebner_basis(generators, (x, y), order="lex")
    direct = sp.groebner(generators, x, y, order="lex")
    assert _exprs(accelerated) == _exprs(direct)


def test_modular_groebner_discards_unlucky_leading_coefficient_prime():
    x, y = sp.symbols("x y")
    generators = (32003 * x + y, x**2 + y**2 - 1)
    result = modular_groebner_basis_qq(generators, (x, y), order="grevlex", max_primes=8)
    assert result is not None
    assert 32003 not in result.certificate.primes
    direct = sp.groebner(generators, x, y, order="grevlex", domain=sp.QQ)
    assert result.basis == _exprs(direct)
