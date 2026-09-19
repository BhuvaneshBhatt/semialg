import pytest
import sympy as sp

from semialg.algebraic import (
    certified_sign_stable_root_neighborhood,
    certify_polynomial_root_interval,
    isolate_real_roots,
    rational_between_algebraic_reals,
    verify_polynomial_root_interval_certificate,
)


def test_descartes_certifies_root_free_and_unique_intervals_without_sturm():
    x = sp.Symbol("x")
    poly = sp.Poly((x - 1) * (x + 2), x, domain=sp.QQ)

    empty = certify_polynomial_root_interval(poly, 2, 3)
    unique = certify_polynomial_root_interval(poly, 0, 2)

    assert empty.root_count == 0
    assert empty.method == "descartes"
    assert unique.root_count == 1
    assert unique.method == "descartes"
    assert verify_polynomial_root_interval_certificate(unique)


def test_ambiguous_descartes_interval_falls_back_to_exact_sturm_count():
    x = sp.Symbol("x")
    poly = sp.Poly((x + 3) * (x + 1) * (x - 1) * (x - 2), x, domain=sp.QQ)
    cert = certify_polynomial_root_interval(poly, -4, 4)
    assert cert.root_count == 4
    assert cert.descartes_variations is not None
    assert cert.descartes_variations > 1
    assert cert.method == "descartes+sturm"
    assert verify_polynomial_root_interval_certificate(cert)


def test_interval_endpoint_semantics_are_exact():
    x = sp.Symbol("x")
    poly = sp.Poly(x * (x - 1), x, domain=sp.QQ)
    closed = certify_polynomial_root_interval(poly, 0, 1)
    open_interval = certify_polynomial_root_interval(
        poly, 0, 1, include_left=False, include_right=False
    )
    assert closed.root_count == 2
    assert closed.left_is_root and closed.right_is_root
    assert open_interval.root_count == 0


def test_rational_separator_between_certified_algebraic_roots():
    x = sp.Symbol("x")
    roots = isolate_real_roots(sp.Poly(x**2 - 2, x, domain=sp.QQ))
    separator = rational_between_algebraic_reals(roots[0], roots[1])
    assert separator.is_Rational
    assert sp.N(roots[0].as_expr()) < separator < sp.N(roots[1].as_expr())


def test_sign_stable_neighborhood_preserves_companion_signs():
    x = sp.Symbol("x")
    positive_root = isolate_real_roots(sp.Poly(x**2 - 2, x, domain=sp.QQ))[1]
    result = certified_sign_stable_root_neighborhood(
        positive_root,
        (sp.Poly(x + 1, x, domain=sp.QQ), sp.Poly(3 - x, x, domain=sp.QQ)),
    )
    assert result.signs == (1, 1)
    assert result.interval.left < result.interval.right
    for companion in result.companion_polynomials:
        cert = certify_polynomial_root_interval(
            companion, result.interval.left, result.interval.right
        )
        assert cert.root_count == 0


def test_sign_stable_neighborhood_rejects_companion_zero_at_root():
    x = sp.Symbol("x")
    root = isolate_real_roots(sp.Poly(x**2 - 2, x, domain=sp.QQ))[1]
    with pytest.raises(ValueError, match="vanishes"):
        certified_sign_stable_root_neighborhood(root, (x**2 - 2,))


def test_certified_root_interval_api_rejects_transcendental_inputs():
    x = sp.Symbol("x")
    with pytest.raises((TypeError, sp.PolynomialError)):
        certify_polynomial_root_interval(sp.sin(x) - x / 2, -1, 1, var=x)


def test_existing_transcendental_solver_subpackage_remains_shipped():
    import importlib.util

    assert importlib.util.find_spec("semialg.solve.transcendental") is not None
