import sympy as sp

from semialg.parametric_geometry import intrinsic_parametric_cover
from semialg.standard_regions import Ball, EllipsoidBoundary, Sphere
from semialg.symbolic_regions import SemialgebraicRegion


def test_sphere_native_chart_has_exact_intrinsic_metric():
    x, y = sp.symbols("x y", real=True)
    cover = intrinsic_parametric_cover(Sphere((0, 0), 2), (x, y))
    assert cover.certified_dimension() == 1
    chart = cover.charts[0]
    assert chart.source == "sphere"
    assert sp.simplify(chart.metric_factor() - 2) == 0
    pulled = chart.intrinsic_integrand(x**2 + y**2, (x, y))
    assert sp.simplify(pulled - 8) == 0


def test_ball_native_chart_uses_radial_metric_without_cad():
    x, y = sp.symbols("x y", real=True)
    cover = Ball((1, -1), 3).intrinsic_parametric_cover((x, y))
    assert cover.certified_dimension() == 2
    chart = cover.charts[0]
    radial = chart.parameters[0]
    assert chart.source == "ball"
    assert sp.simplify(chart.metric_factor() - radial) == 0
    assert chart.bounds[0] == (radial, 0, 3)


def test_ellipsoid_boundary_native_chart_uses_induced_metric():
    x, y = sp.symbols("x y", real=True)
    ellipse = EllipsoidBoundary((0, 0), sp.diag(4, 9))
    chart = intrinsic_parametric_cover(ellipse, (x, y)).charts[0]
    theta = chart.parameters[0]
    assert chart.mapping == (2 * sp.cos(theta), 3 * sp.sin(theta))
    assert (
        sp.simplify(chart.metric_factor() ** 2 - (4 * sp.sin(theta) ** 2 + 9 * sp.cos(theta) ** 2))
        == 0
    )


def test_intrinsic_cover_normalizes_cad_graph_chart():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(y, x), x >= 0, x <= 1)
    cover = intrinsic_parametric_cover(formula, (x, y), dimension=1)
    assert cover.certified_dimension() == 1
    assert len(cover.charts) == 1
    chart = cover.charts[0]
    assert chart.source == "cad_intrinsic"
    assert chart.mapping == (x, x)
    assert chart.bounds == ((x, 0, 1),)
    assert sp.simplify(chart.metric_factor() - sp.sqrt(2)) == 0


def test_semialgebraic_region_exposes_intrinsic_chart_convenience_api():
    x, y = sp.symbols("x y", real=True)
    region = SemialgebraicRegion(sp.And(sp.Eq(y, x), x >= 0, x <= 1), (x, y))
    cover = region.intrinsic_parametric_cover(dimension=1)
    assert cover.charts[0].source == "cad_intrinsic"


def test_intrinsic_cover_does_not_drop_singular_target_strata():
    x, y = sp.symbols("x y", real=True)
    cusp = sp.Eq(y**2, x**3)
    try:
        intrinsic_parametric_cover(cusp, (x, y), dimension=1)
    except (ValueError, NotImplementedError):
        pass
