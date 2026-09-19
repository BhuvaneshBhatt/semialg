import sympy as sp

from semialg.algebraic import (
    ModularResultantCertificate,
    ModularSubresultantCertificate,
    modular_resultant_qq,
    modular_subresultants_qq,
    subresultant_prs,
    verify_modular_resultant_certificate,
    verify_modular_subresultant_certificate,
)
from semialg.cad_algorithms.polynomial_utils import resultant_expression


def test_modular_resultant_reconstructs_parameter_polynomial():
    a, b, x = sp.symbols("a b x")
    f = x**2 + a * x + b
    g = x + a
    result = modular_resultant_qq(f, g, x)
    assert result is not None
    assert result.resultant == b
    assert len(result.certificate.primes) >= 2
    assert verify_modular_resultant_certificate(result.certificate)


def test_modular_resultant_handles_rational_coefficients():
    a, x = sp.symbols("a x")
    f = sp.Rational(2, 3) * x**2 + a
    g = x + sp.Rational(5, 7)
    result = modular_resultant_qq(f, g, x)
    assert result is not None
    assert sp.expand(result.resultant - (a + sp.Rational(50, 147))) == 0
    assert verify_modular_resultant_certificate(result.certificate)


def test_modular_resultant_discards_unlucky_degree_drop_prime():
    a, b, x = sp.symbols("a b x")
    f = 32003 * x**2 + a * x + 1
    g = x + b
    result = modular_resultant_qq(f, g, x, max_primes=8)
    assert result is not None
    assert 32003 not in result.certificate.primes
    assert sp.expand(result.resultant - (32003 * b**2 - a * b + 1)) == 0


def test_resultant_certificate_rejects_tampering():
    a, x = sp.symbols("a x")
    result = modular_resultant_qq(x**2 + a, x + 1, x)
    assert result is not None
    cert = result.certificate
    bad = ModularResultantCertificate(
        cert.variable,
        cert.parameters,
        cert.first,
        cert.second,
        sp.expand(cert.resultant + 1),
        cert.parameter_degree_bounds,
        cert.verification_grid,
        cert.primes,
        cert.modulus,
    )
    assert not verify_modular_resultant_certificate(bad)


def test_modular_subresultants_reconstruct_exact_native_normalization():
    a, b, x = sp.symbols("a b x")
    f = x**3 + a * x + b
    g = x**2 + a
    result = modular_subresultants_qq(f, g, x)
    assert result is not None
    direct = tuple(
        sp.expand(poly.as_expr() if isinstance(poly, sp.Poly) else poly)
        for poly in sp.subresultants(sp.Poly(f, x), sp.Poly(g, x))
    )
    assert result.polynomials == direct
    assert verify_modular_subresultant_certificate(result.certificate)


def test_subresultant_certificate_rejects_tampering():
    a, x = sp.symbols("a x")
    result = modular_subresultants_qq(x**2 + a * x + 1, x + a, x)
    assert result is not None
    cert = result.certificate
    polys = list(cert.polynomials)
    polys[-1] = sp.expand(polys[-1] + 1)
    bad = ModularSubresultantCertificate(
        cert.variable,
        cert.parameters,
        cert.first,
        cert.second,
        tuple(polys),
        cert.resultant_certificate,
        cert.primes,
        cert.modulus,
    )
    assert not verify_modular_subresultant_certificate(bad)


def test_existing_subresultant_helper_uses_modular_acceleration():
    a, x = sp.symbols("a x")
    result = subresultant_prs(x**2 + a * x + 1, x + a, x)
    assert result.source.startswith("multimodular_subresultants")
    assert result.resultant == 1


def test_cad_resultant_helper_matches_direct_exact_resultant():
    a, x = sp.symbols("a x")
    left = sp.Poly(x**2 + a * x + 1, x)
    right = sp.Poly(x + a, x)
    assert resultant_expression(left, right, x) == sp.resultant(left.as_expr(), right.as_expr(), x)


def test_modular_resultant_certifies_zero_for_common_factor():
    a, x = sp.symbols("a x")
    common = x + a
    result = modular_resultant_qq(common * (x + 1), common * (x + 2), x)
    assert result is not None
    assert result.resultant == 0
    assert verify_modular_resultant_certificate(result.certificate)


def test_modular_resultant_declines_excessive_verification_grid():
    a, b, x = sp.symbols("a b x")
    f = x**4 + a**5 * x + b**5
    g = x**4 + a**5 + b**5 * x
    assert modular_resultant_qq(f, g, x, max_verification_points=4) is None
