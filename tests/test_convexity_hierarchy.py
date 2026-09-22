from __future__ import annotations

import sympy as sp

from semialg import convexity_certificate, is_convex
from semialg.convexity import (
    ConvexityCertificate,
    PolynomialConvexityCertificate,
    QuadraticConvexityCertificate,
    polynomial_convexity_certificate,
    quadratic_convexity_certificate,
)


def test_1d_classifier_rejects_disconnected_relation():
    x = sp.symbols("x", real=True)
    cert = convexity_certificate(x**2 >= 1, (x,))

    assert cert.outcome is False
    assert cert.method == "one-dimensional-ordered-cad-interval"
    assert cert.details["component_count"] == 2


def test_complete_one_dimensional_classifier_accepts_open_interval():
    x = sp.symbols("x", real=True)
    cert = convexity_certificate(sp.And(x > -2, x < 3), (x,))

    assert cert.outcome is True
    assert cert.method == "one-dimensional-ordered-cad-interval"


def test_empty_set_has_trivial_convexity_certificate():
    x, y = sp.symbols("x y", real=True)
    cert = convexity_certificate(sp.false, (x, y))

    assert cert == ConvexityCertificate(True, "empty-set", sp.false, (x, y))


def test_nonlinear_singleton_is_recognized_before_general_fallback():
    x, y = sp.symbols("x y", real=True)
    cert = convexity_certificate(sp.Eq(x**2 + y**2, 0), (x, y))

    assert cert.outcome is True
    assert cert.method == "singleton"
    assert cert.witness == {x: 0, y: 0}


def test_quadratic_classifier_certifies_ellipsoid():
    x, y = sp.symbols("x y", real=True)
    cert = quadratic_convexity_certificate(3 * x**2 + 2 * x * y + 4 * y**2 <= 1, (x, y))

    assert isinstance(cert, QuadraticConvexityCertificate)
    assert cert.outcome is True
    assert cert.method == "quadratic-intersection"


def test_quadratic_classifier_declines_indefinite_sublevel():
    x, y = sp.symbols("x y", real=True)
    cert = quadratic_convexity_certificate(x**2 - y**2 <= 1, (x, y))

    assert cert.outcome is None
    assert cert.method == "indefinite-quadratic"


def test_polynomial_certificate_infers_variables_and_proves_global_convexity():
    x, y = sp.symbols("x y", real=True)
    cert = polynomial_convexity_certificate(x**4 + y**4)

    assert isinstance(cert, PolynomialConvexityCertificate)
    assert cert.variables == (x, y)
    assert cert.outcome is True


def test_domain_relative_hessian_can_succeed_when_global_hessian_fails():
    x, y = sp.symbols("x y", real=True)
    expr = x**4 - x**2 + y**2

    global_cert = polynomial_convexity_certificate(expr, (x, y))
    domain_cert = polynomial_convexity_certificate(expr, (x, y), domain=x >= 1)
    region_cert = convexity_certificate(sp.And(x >= 1, expr <= 10), (x, y))

    assert global_cert.outcome is False
    assert domain_cert.outcome is True
    assert region_cert.outcome is True
    assert region_cert.method == "domain-relative-hessian"


def test_topology_stage_rejects_disconnected_factorized_variety():
    x, y = sp.symbols("x y", real=True)
    circles = sp.Eq((x**2 + y**2 - 1) * (x**2 + y**2 - 4), 0)
    cert = convexity_certificate(circles, (x, y))

    assert cert.outcome is False
    assert cert.method in {"disconnected-factorized-variety", "disconnected-cad"}


def test_exact_segment_search_rejects_connected_l_shape():
    x, y = sp.symbols("x y", real=True)
    horizontal = sp.And(x >= 0, x <= 2, y >= 0, y <= 1)
    vertical = sp.And(x >= 0, x <= 1, y >= 0, y <= 2)
    cert = convexity_certificate(sp.Or(horizontal, vertical), (x, y))

    assert cert.outcome is False
    assert cert.method == "exact-segment-counterexample"
    left = cert.details["left"]
    right = cert.details["right"]
    midpoint = cert.witness
    region = sp.Or(horizontal, vertical)
    assert bool(region.subs(left))
    assert bool(region.subs(right))
    assert not bool(region.subs(midpoint))


def test_is_convex_delegates_to_certificate_hierarchy():
    x, y = sp.symbols("x y", real=True)
    region = x**2 + y**2 <= 1

    cert = convexity_certificate(region, (x, y))
    assert cert.outcome is True
    assert cert.method == "quadratic-intersection"
    assert is_convex(region, (x, y)) is True


def test_polynomial_hessian_sign_is_vacuously_true_on_empty_domain():
    x = sp.symbols("x", real=True)
    cert = polynomial_convexity_certificate(-(x**2), (x,), domain=x**2 < 0)
    assert cert.outcome is True


def test_affine_equality_region_uses_affine_stage_without_singleton_qe():
    x, y = sp.symbols("x y", real=True)
    cert = convexity_certificate(sp.Eq(x, 0), (x, y))
    assert cert.outcome is True
    assert cert.method == "affine-polyhedron"


def test_high_dimensional_quadratic_uses_ldlt_inertia_without_principal_minors():
    xs = sp.symbols("x0:8", real=True)
    expr = sum((i + 1) * xs[i] ** 2 for i in range(len(xs)))
    cert = polynomial_convexity_certificate(expr, xs)

    assert cert.outcome is True
    assert cert.method == "constant-hessian-ldlt-inertia"
    assert cert.principal_minors == ()
    assert cert.hessian_inertia == (8, 0, 0)


def test_constant_hessian_inertia_detects_indefinite_quadratic():
    x, y, z = sp.symbols("x y z", real=True)
    cert = polynomial_convexity_certificate(x**2 + y**2 - z**2, (x, y, z))

    assert cert.outcome is False
    assert cert.method == "constant-hessian-ldlt-inertia"
    assert cert.hessian_inertia == (2, 1, 0)


def test_convexity_cad_analysis_reuses_one_cylindrical_solution():
    import semialg.convexity as convexity_module

    x, y = sp.symbols("x y", real=True)
    circle = sp.Eq(x**2 + y**2, 1)
    calls = 0
    original = convexity_module.extract_cylindrical_solution

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    analysis = convexity_module._ConvexityCADAnalysis(circle, (x, y), extractor=counted)
    first = analysis.cylindrical_solution()
    second = analysis.cylindrical_solution()
    analysis.connectivity_graph()

    assert first is second
    assert calls == 1


def test_selected_cad_samples_are_lightweight_exact_representatives():
    from semialg.cad_algorithms import extract_selected_cad_samples

    x, y = sp.symbols("x y", real=True)
    region = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    samples = extract_selected_cad_samples(region, (x, y), limit=3)

    assert 1 <= len(samples) <= 3
    assert all(bool(region.subs(sample.point)) for sample in samples)
    assert all(len(sample.index) == 2 for sample in samples)


def test_polynomial_convexity_certificate_cache_reuses_exact_result():
    from semialg.convexity import _cached_polynomial_convexity_certificate

    x, y = sp.symbols("x y", real=True)
    expr = x**4 - x**2 + y**2
    _cached_polynomial_convexity_certificate.cache_clear()
    first = polynomial_convexity_certificate(expr, (x, y), domain=x >= 1)
    before = _cached_polynomial_convexity_certificate.cache_info()
    second = polynomial_convexity_certificate(expr, (x, y), domain=x >= 1)
    after = _cached_polynomial_convexity_certificate.cache_info()

    assert first == second
    assert after.hits == before.hits + 1


def test_factorized_definite_quadrics_use_early_disconnected_rejection():
    x, y = sp.symbols("x y", real=True)
    circles = sp.Eq((x**2 + y**2 - 1) * (x**2 + y**2 - 4), 0)

    cert = convexity_certificate(circles, (x, y))

    assert cert.outcome is False
    assert cert.method == "disconnected-factorized-variety"
