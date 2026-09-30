import sympy as sp

from semialg import (
    blowup_charts,
    contact_order,
    curve_selection,
    local_bound,
    local_components,
    local_geometry,
    local_germ,
    local_image,
    local_preimage,
    local_range,
    local_sign_strata,
    parameter_strata,
    path_independent,
    vanishing_order,
)


def test_cusp_geometry_for_asymptotic_consumer():
    x, y = sp.symbols("x y", real=True)
    region = sp.And(sp.Eq(y, x**2), x > 0)
    geometry = local_geometry(region, (0, 0), (x, y))
    curve = curve_selection(region, (0, 0), (x, y))
    assert geometry.dimension == 1
    assert (vanishing_order(x, curve), vanishing_order(y, curve)) == (1, 2)


def test_branch_and_sign_information_for_limit_consumer():
    x = sp.Symbol("x", real=True)
    assert len(local_components(sp.Ne(x, 0), (0,), (x,))) == 2
    assert {piece.signs for piece in local_sign_strata((x,), sp.Ne(x, 0), (0,), (x,))} == {
        (-1,),
        (1,),
    }


def test_weighted_chart_exposes_relative_scale():
    x, y = sp.symbols("x y", real=True)
    chart = blowup_charts((x, y), weights=(1, 2))[0]
    assert sp.simplify(chart.mapping[1] / chart.mapping[0] ** 2) == chart.directions[0]


def test_local_certification_surface_for_asymptotic_consumer():
    x, y, a = sp.symbols("x y a", real=True)
    germ = local_germ(sp.Eq(y, x**2), (0, 0), (x, y))
    assert path_independent(y - x**2, germ, 0)
    assert local_range(x, sp.And(x >= -1, x <= 1), (0,), (x,), radius=1).upper == 1
    assert local_bound(x, sp.true, (0,), (x,), bound=1, return_certificate=True).certified
    image = local_image((x,), local_germ(x >= 0, (0,), (x,)), image_variables=(y,))
    assert local_preimage(image, (x,), (x,)).point[x] == 0
    assert len(parameter_strata(sp.Eq(x**2, a), (x,), (a,))) == 2
    curve = curve_selection(sp.And(sp.Eq(y, x**2), x > 0), (0, 0), (x, y))
    assert contact_order(y, 0, curve) == 2


def test_local_geometry_root_api_regression_contracts():
    """Exercise the local/asymptotic root APIs through one compact exact example."""
    x, y = sp.symbols("x y", real=True)
    region = sp.And(sp.Eq(y, x**2), x > 0)
    geometry = local_geometry(region, (0, 0), (x, y))
    curve = curve_selection(region, (0, 0), (x, y))
    assert geometry.dimension == 1
    assert vanishing_order(y, curve) == 2
    assert contact_order(y, 0, curve) == 2
    assert len(blowup_charts((x, y), weights=(1, 2))) == 2
    germ = local_germ(region, (0, 0), (x, y))
    assert path_independent(y - x**2, germ, 0)
    assert local_range(x, sp.And(x >= -1, x <= 1), (0,), (x,), radius=1).upper == 1
    image = local_image((x,), local_germ(x >= 0, (0,), (x,)), image_variables=(y,))
    assert image.point[y] == 0
