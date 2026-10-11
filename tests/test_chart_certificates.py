"""Geometry certificates never infer integrability or accept unproved null fibers."""

from dataclasses import replace

import pytest
import sympy as sp

from semialg.chart_certificates import (
    Chart,
    GeometryBudget,
    GeometryProofError,
    certify_chart,
    certify_cover,
    certify_linear_metric,
    certify_null_set,
    chart_cache_info,
    clear_chart_cache,
    polar_chart,
)

x, y, t, u, a = sp.symbols("x y t u a", real=True)


@pytest.mark.parametrize("condition", [a > 0, a < 0, sp.Ne(a, 0)])
def test_null_polynomial_nonzero_on_every_parameter_fiber(condition):
    c = certify_null_set(sp.Eq(a * x, 0), (x, y), assumptions=condition)
    assert c.zero_polynomials == (a * x,)


@pytest.mark.parametrize("region", [sp.Eq(a * x, 0), y < 1, sp.true])
def test_positive_volume_or_exceptional_parameter_null_claim_rejected(region):
    with pytest.raises(GeometryProofError):
        certify_null_set(region, (x, y))


def test_rank_loss_parameter_cannot_be_certified_generically():
    c = Chart((t, u), (a * t, u), (x / a, y), ((u, 0, 1), (t, 0, 1)))
    assert certify_chart(c, x / a, (x, y), assumptions=a > 0).determinant == a
    with pytest.raises(GeometryProofError):
        certify_chart(c, x / a, (x, y))


def test_cover_missing_component_overlap_and_forged_certificate():
    c = Chart((t, u), (t, u), (x, y), ((u, 0, 1), (t, 0, 1)))
    proof = certify_chart(c, x, (x, y))
    region = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    assert certify_cover(region, (proof,)).source_null.region is sp.false
    with pytest.raises(GeometryProofError):
        certify_cover(region, (proof, proof))
    with pytest.raises(GeometryProofError):
        certify_cover(region, (replace(proof, determinant=3),))
    with pytest.raises(GeometryProofError):
        certify_cover(sp.And(x >= 0, x <= 2, y >= 0, y <= 1), (proof,))


def test_cache_key_includes_assumptions_and_budget_and_is_bounded():
    clear_chart_cache()
    c = Chart((t, u), (t, u), (x, y), ((u, 0, 1), (t, 0, 1)))
    p = certify_chart(c, x, (x, y))
    assert certify_chart(c, x, (x, y)) is p
    assert chart_cache_info()[1].hits == 1
    with pytest.raises(GeometryProofError):
        certify_chart(c, x, (x, y), budget=GeometryBudget(max_operations=0))
    assert chart_cache_info()[1].maxsize == 128
    clear_chart_cache()
    assert all(i.currsize == 0 for i in chart_cache_info())


@pytest.mark.parametrize("squared", [False, True])
def test_polar_null_seam_and_coverage(squared):
    q = x * x + y * y
    c = polar_chart((x, y), (t, u), (0, 1), squared=squared)
    proof = certify_chart(c, q if squared else sp.sqrt(q), (x, y))
    cover = certify_cover(q <= 1, (proof,))
    assert cover.source_null.zero_polynomials == (y,)
    assert proof.map_jacobian_squared == 4 * t if squared else proof.map_jacobian_squared == 1


def test_null_regular_chart_rejects_positive_volume_exclusions():
    c = Chart((t, u), (t * u, u), (x / y, y), ((u, 0, 1), (t, 0, 1)), sp.Eq(u, 0), sp.Eq(y, 0))
    p = certify_chart(c, x / y, (x, y))
    assert p.determinant == u and p.chart_null.zero_polynomials == (u,)
    with pytest.raises(GeometryProofError):
        certify_chart(replace(c, chart_exclusions=u < sp.Rational(1, 2)), x / y, (x, y))


def test_algebraic_positive_branch_and_wrong_branch():
    c = Chart(
        (t, u),
        (t, sp.sqrt(u)),
        (x, y * y),
        ((u, 0, 1), (t, 0, 1)),
        sp.Eq(u, 0),
        sp.Eq(y, 0),
        y >= 0,
        "algebraic",
    )
    p = certify_chart(c, x, (x, y))
    assert sp.simplify(p.determinant - 1 / (2 * sp.sqrt(u))) == 0
    with pytest.raises(GeometryProofError):
        certify_chart(replace(c, source_domain=sp.true), x, (x, y))


@pytest.mark.parametrize(
    "matrix,expected", [(sp.eye(2), 1), (sp.Matrix([[1, 0], [0, 1], [2, 3]]), 14)]
)
def test_shared_linear_metric_for_ambient_and_intrinsic_charts(matrix, expected):
    assert certify_linear_metric(matrix).gram_determinant == expected


def test_linear_metric_parameter_rank_loss_and_zero_map():
    assert (
        certify_linear_metric([[a, 0], [0, 1]], assumptions=sp.Ne(a, 0)).gram_determinant == a * a
    )
    for matrix in [[[a, 0], [0, 1]], [[0, 0], [0, 1]]]:
        with pytest.raises(GeometryProofError):
            certify_linear_metric(matrix)


def test_geometry_guide_examples_execute():
    import re
    from pathlib import Path

    namespace = {}
    guide = (Path(__file__).parents[1] / "docs/chart_certificates.md").read_text()
    for code in re.findall(r"```python\n(.*?)```", guide, re.S):
        exec(compile(code, "chart_certificates.md", "exec"), namespace)


def test_cached_geometry_ignores_unrecorded_global_parameter_and_branch_facts():
    offset = sp.Symbol("offset")
    shifted = Chart((t, u), (t + offset, u), (x - offset, y), ((u, 0, 1), (t, 0, 1)))
    clear_chart_cache()
    with sp.assuming(sp.Q.real(offset)):
        with pytest.raises(GeometryProofError):
            certify_chart(shifted, x - offset, (x, y))
    proved = certify_chart(shifted, x - offset, (x, y), assumptions=sp.Q.real(offset))
    assert proved.assumptions == (sp.Q.real(offset),)
    with pytest.raises(GeometryProofError):
        certify_chart(shifted, x - offset, (x, y))
    root = Chart(
        (t, u),
        (t, sp.sqrt(u)),
        (x, y * y),
        ((u, 0, 1), (t, 0, 1)),
        sp.Eq(u, 0),
        sp.Eq(y, 0),
        kind="algebraic",
    )
    with sp.assuming(sp.Q.positive(y)):
        with pytest.raises(GeometryProofError):
            certify_chart(root, x, (x, y))


def test_explicit_Q_parameter_sign_and_complex_nonzero_distinction():
    c = Chart((t, u), (a * t, u), (x / a, y), ((u, 0, 1), (t, 0, 1)))
    assert certify_chart(c, x / a, (x, y), assumptions=sp.Q.positive(a)).determinant == a
    offset = sp.Symbol("complex_offset")
    shifted = Chart((t, u), (t + offset, u), (x - offset, y), ((u, 0, 1), (t, 0, 1)))
    with pytest.raises(GeometryProofError):
        certify_chart(shifted, x - offset, (x, y), assumptions=sp.Ne(offset, 0))


def test_chart_growth_budget_precedes_polynomial_expansion(monkeypatch):
    c = Chart(
        (t, u), (t + (u + 1) ** 1000000, u), (x - (y + 1) ** 1000000, y), ((u, 0, 1), (t, 0, 1))
    )
    monkeypatch.setattr(sp, "cancel", lambda *a, **k: pytest.fail("unbounded expansion"))
    with pytest.raises(GeometryProofError, match="before expansion"):
        certify_chart(c, x - (y + 1) ** 1000000, (x, y))


def test_null_polynomial_requires_finite_real_coefficients():
    with pytest.raises(GeometryProofError):
        certify_null_set(sp.Eq(sp.oo * x, 0, evaluate=False), (x, y))


def test_vector_map_gram_identity_and_rank_loss():
    from semialg.chart_certificates import Chart, GeometryProofError, certify_chart

    x, y, z, t, s, u = sp.symbols("x y z t s u", real=True)
    chart = Chart(
        (t, s, u), (t / 2, s / 3, u), (2 * x, 3 * y, z), ((u, 0, 1), (s, 0, 3), (t, 0, 2))
    )
    c = certify_chart(chart, (2 * x, 3 * y), (x, y, z))
    assert (
        c.map_jacobian_squared == 36
        and c.fiber_metric_squared == 1
        and c.determinant == sp.Rational(1, 6)
    )
    with pytest.raises(GeometryProofError):
        certify_chart(chart, (x, x), (x, y, z))


def test_chart_parameter_events_are_candidates_not_generic_certificates():
    from semialg.chart_certificates import (
        Chart,
        GeometryBudget,
        GeometryProofError,
        chart_parameter_events,
    )

    x, y, t, u, a, b = sp.symbols("x y t u a b", real=True)
    chart = Chart((t, u), (a * t, u), (x / a, y), ((u, 0, b), (t, 0, 1)))
    assert chart_parameter_events(chart, x / a, (x, y), (a, b)) == (a, b)
    with pytest.raises(GeometryProofError):
        chart_parameter_events(chart, x / a, (x, y), (a, b), budget=GeometryBudget(max_variables=1))


def test_determinant_expansion_budget_precedes_dense_gram_expansion(monkeypatch):
    from semialg import chart_certificates as gc

    symbols = sp.symbols("c0:8", real=True)
    dense = sp.Matrix(3, 3, lambda i, j: sum(symbols[(i + j + k) % 8] for k in range(5)))

    def forbidden(*args, **kwargs):
        pytest.fail("over-budget determinant expanded")

    monkeypatch.setattr(type(dense), "det", forbidden)
    with pytest.raises(gc.GeometryProofError, match="expansion budget"):
        gc._determinant(dense, gc.GeometryBudget(max_operations=100))


def test_parameter_endpoint_events_include_coincident_shell_boundaries():
    from semialg.chart_certificates import parameter_geometry_events

    a, b = sp.symbols("a b", real=True)
    assert parameter_geometry_events((a - b, 1 / (a + 1), 1 - a * a), (a, b)) == (
        a - 1,
        a + 1,
        a - b,
    )


def test_coordinate_sign_flags_are_not_closed_domain_proofs():
    from semialg.chart_certificates import Chart, GeometryProofError, certify_chart

    x, t = sp.symbols("scope_x scope_t", real=True)
    y, u = sp.symbols("scope_y scope_u", positive=True)
    chart = Chart((t, u), (t * u, u), (x / y, y), ((u, -1, 1), (t, 0, 1)))
    with pytest.raises(GeometryProofError, match="plain real"):
        certify_chart(chart, x / y, (x, y))
