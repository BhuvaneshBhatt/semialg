import sympy as sp

from semialg import (
    BlowupChart,
    CurveSelectionWitness,
    LocalBoundCertificate,
    LocalGeometry,
    LocalGerm,
    LocalRangeResult,
    LocalSignStratum,
    ParameterStratum,
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


def test_local_geometry_line():
    x, y = sp.symbols("x y", real=True)
    result = local_geometry(sp.Eq(y, 0), (0, 0), (x, y))
    assert isinstance(result, LocalGeometry)
    assert isinstance(result.germ, LocalGerm)
    assert result.dimension == 1
    assert sp.simplify(result.tangent_cone.subs(y, 0)) is sp.true


def test_local_sign_strata_crossing():
    x = sp.Symbol("x", real=True)
    strata = local_sign_strata((x,), sp.true, (0,), (x,))
    assert all(isinstance(piece, LocalSignStratum) for piece in strata)
    assert {piece.signs for piece in strata} == {(-1,), (0,), (1,)}


def test_local_components_punctured_line():
    x = sp.Symbol("x", real=True)
    components = local_components(sp.true, (0,), (x,))
    assert len(components) == 2


def test_curve_selection_and_orders():
    x, y = sp.symbols("x y", real=True)
    curve = curve_selection(sp.And(sp.Eq(y, x**2), x > 0), (0, 0), (x, y))
    assert isinstance(curve, CurveSelectionWitness)
    assert vanishing_order(x, curve) == 1
    assert vanishing_order(y, curve) == 2
    assert contact_order(y, 0, curve) == 2


def test_local_bound_interval():
    x = sp.Symbol("x", real=True)
    result = local_range(x, sp.true, (0,), (x,), radius=2)
    assert isinstance(result, LocalRangeResult)
    assert local_bound(x, sp.true, (0,), (x,), radius=2) == 2


def test_local_bound_finds_radius():
    x = sp.Symbol("x", real=True)
    cert = local_bound(x, sp.true, (0,), (x,), bound=1, return_certificate=True)
    assert isinstance(cert, LocalBoundCertificate)
    assert cert.certified
    assert cert.bound == 1
    assert cert.radius.is_positive is True or cert.radius > 0


def test_parameter_feasibility_strata():
    x, a = sp.symbols("x a", real=True)
    strata = parameter_strata(sp.Eq(x**2, a), (x,), (a,))
    assert len(strata) == 2
    assert all(isinstance(piece, ParameterStratum) for piece in strata)
    assert strata[0].value is True
    assert strata[1].value is False


def test_weighted_blowup_chart():
    x, y = sp.symbols("x y", real=True)
    charts = blowup_charts((x, y), weights=(1, 2))
    assert len(charts) == 2
    assert all(isinstance(chart, BlowupChart) for chart in charts)
    assert charts[0].mapping[0] == charts[0].radial
    assert charts[0].mapping[1].has(charts[0].radial ** 2)


def test_local_image_preimage():
    x, y = sp.symbols("x y", real=True)
    germ = local_germ(x >= 0, (0,), (x,))
    image = local_image((x**2,), germ, image_variables=(y,))
    assert image.point[y] == 0
    preimage = local_preimage(image, (x,), (x,))
    assert preimage.point[x] == 0


def test_path_independence():
    x, y = sp.symbols("x y", real=True)
    germ = local_germ(sp.Eq(y, x**2), (0, 0), (x, y))
    assert path_independent(y - x**2, germ, 0)
    assert not path_independent(x, germ, 0)


def test_curve_selection_theorem_fallback():
    x, y = sp.symbols("x y", real=True)
    result = curve_selection(sp.And(sp.Eq(y, x + x**2), x > 0), (0, 0), (x, y), max_order=2)
    assert result.certified
    assert not result.explicit
    assert result.method == "curve-selection-theorem"
    assert result.existence_formula is not None


def test_local_bound_finds_bound_and_radius():
    x = sp.Symbol("x", real=True)
    cert = local_bound(x, sp.true, (0,), (x,), return_certificate=True)
    assert cert.certified
    assert cert.radius > 0
    assert cert.bound >= 0


def test_local_geometry_result_types_construct():
    x, t = sp.symbols("x t", real=True)
    germ = LocalGerm(sp.true, {x: 0}, (x,))
    assert LocalGeometry(germ, 1, sp.true).dimension == 1
    assert LocalSignStratum(x > 0, (1,)).signs == (1,)
    assert CurveSelectionWitness(t, (t,), t > 0, {x: 0}, (x,), x > 0).explicit
    assert LocalRangeResult(-1, 1, True, True, sp.Abs(x) <= 1).upper == 1
    assert LocalBoundCertificate(1, 1, sp.Abs(x) < 1, sp.true).certified
    assert ParameterStratum(x > 0, True).value is True
    assert BlowupChart(t, (), (t,), (1,), 0).pivot == 0
