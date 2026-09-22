import pytest
import sympy as sp

from semialg import (
    Ball,
    BooleanRegion,
    Box,
    Capsule,
    Cylinder,
    Interval,
    Parallelepiped,
    Parallelogram,
    ParametricRegion,
    Polygon,
    Prism,
    Pyramid,
    Simplex,
    Sphere,
    SphericalShell,
    Stadium,
)
from semialg.parametric_integration import (
    integrate_over_parametric_region,
    reduce_parametric_integral,
)
from semialg.standard_region_integrate import integrate_over_standard_region


def test_parametric_integral_preserves_ambient_symbol_identity():
    x = sp.Symbol("x")
    u = sp.Symbol("u")
    region = ParametricRegion((u,), ((u, 0, 1),), (u,))
    reduced, _, _ = reduce_parametric_integral(x, ["x"], region)
    assert reduced == u
    assert integrate_over_parametric_region(x, ["x"], region) == sp.Rational(1, 2)


def test_parametric_region_requires_exact_parameter_limits():
    u, v = sp.symbols("u v")
    with pytest.raises(ValueError, match="not a declared parameter"):
        ParametricRegion((u,), ((v, 0, 1),), (u,))
    with pytest.raises(ValueError, match="not a declared parameter"):
        ParametricRegion((u,), (("v", 0, 1),), (u,))
    with pytest.raises(ValueError, match="missing integration limits"):
        ParametricRegion((u,), (), (u,))
    with pytest.raises(ValueError, match="duplicate integration limit"):
        ParametricRegion((u,), ((u, 0, 1), (u, 0, 1)), (u,))


def test_parametric_region_requires_positive_multiplicity():
    u = sp.Symbol("u", real=True)
    with pytest.raises(ValueError, match="positive"):
        ParametricRegion((u,), ((u, 0, 1),), (u,), multiplicity=0)
    with pytest.raises(ValueError, match="positive"):
        ParametricRegion((u,), ((u, 0, 1),), (u,), multiplicity=-1)
    m = sp.Symbol("m", real=True)
    with pytest.raises(ValueError, match="provably positive"):
        ParametricRegion((u,), ((u, 0, 1),), (u,), multiplicity=m)


def test_interval_and_box_reject_reversed_bounds():
    x = sp.Symbol("x", real=True)
    with pytest.raises(ValueError, match="lower endpoint"):
        Interval(2, 1)
    with pytest.raises(ValueError, match="lower endpoint"):
        Box(((0, 1), (2, 1)))
    assert integrate_over_standard_region(1, Interval(1, 2), [x]) == 1


def test_radial_regions_reject_invalid_radii():
    constructors = (
        lambda: Ball((0, 0), -1),
        lambda: Sphere((0, 0, 0), -1),
        lambda: Cylinder((0, 0, 0), (0, 0, 1), -1),
        lambda: Stadium((0, 0), (1, 0), -1),
        lambda: Capsule((0, 0, 0), (0, 0, 1), -1),
    )
    for make_region in constructors:
        with pytest.raises(ValueError, match="nonnegative"):
            make_region()
    with pytest.raises(ValueError, match="shell radii"):
        SphericalShell((0, 0, 0), (2, 1))
    with pytest.raises(ValueError, match="nonnegative"):
        SphericalShell((0, 0, 0), (-1, 2))


def test_standard_regions_validate_ambient_dimensions():
    with pytest.raises(ValueError, match="same dimension"):
        Simplex(((0, 0), (1, 0, 0)))
    with pytest.raises(ValueError, match="origin dimension"):
        Parallelogram((0, 0), ((1, 0), (0, 1, 0)))
    with pytest.raises(ValueError, match="origin dimension"):
        Parallelepiped((0, 0, 0), ((1, 0, 0), (0, 1)))
    base = Polygon(((0, 0), (1, 0), (0, 1)))
    with pytest.raises(ValueError, match="base ambient dimension"):
        Prism(base, (0, 0, 1))
    with pytest.raises(ValueError, match="base ambient dimension"):
        Pyramid(base, (0, 0, 1))
    with pytest.raises(ValueError, match="same dimension"):
        Cylinder((0, 0), (0, 0, 1), 1)


def test_interval_measure_identity_for_small_rational_intervals():
    x = sp.Symbol("x", real=True)
    endpoints = tuple(sp.Rational(i, 2) for i in range(-2, 3))
    for lower in endpoints:
        for upper in endpoints:
            if lower <= upper:
                region = Interval(lower, upper)
                expected = sp.Integer(1) if lower == upper else upper - lower
                assert integrate_over_standard_region(1, region, [x]) == expected


def test_interval_boolean_measure_identities():
    x = sp.Symbol("x", real=True)
    cases = (
        (Interval(0, 2), Interval(1, 3)),
        (Interval(0, 1), Interval(2, 3)),
        (Interval(-1, 3), Interval(0, 1)),
    )
    for left, right in cases:
        m_left = integrate_over_standard_region(1, left, [x])
        m_right = integrate_over_standard_region(1, right, [x])
        m_inter = integrate_over_standard_region(1, BooleanRegion.intersection(left, right), [x])
        m_union = integrate_over_standard_region(1, BooleanRegion.union(left, right), [x])
        assert sp.simplify(m_union + m_inter - m_left - m_right) == 0
        m_diff = integrate_over_standard_region(1, BooleanRegion.difference(left, right), [x])
        assert sp.simplify(m_diff - (m_left - m_inter)) == 0


def test_region_difference_of_disjoint_intervals_preserves_left_region():
    x = sp.Symbol("x", real=True)
    left = Interval(0, 1)
    right = Interval(2, 3)
    assert integrate_over_standard_region(1, BooleanRegion.difference(left, right), [x]) == 1


def test_degenerate_standard_region_dimensions_are_not_overstated():
    assert Box(((0, 0), (0, 1))).dimension() == 1
    assert Simplex(((0, 0), (1, 0), (2, 0))).dimension() == 1
    assert Parallelogram((0, 0), ((1, 0), (2, 0))).dimension() == 1
    assert Ball((0, 0), 0).dimension() == 0
    assert Sphere((0, 0, 0), 0).dimension() == 0
    assert SphericalShell((0, 0, 0), (1, 1)).dimension() == 2
    assert Cylinder((0, 0, 0), (0, 0, 1), 0).dimension() == 1
    assert Capsule((0, 0, 0), (0, 0, 1), 0).dimension() == 1


def test_stadium_rejects_nonplanar_endpoints_but_capsule_allows_them():
    with pytest.raises(ValueError, match="two-dimensional"):
        Stadium((0, 0, 0), (1, 0, 0), 1)
    assert Capsule((0, 0, 0), (1, 0, 0), 1).ambient_dimension() == 3


def test_polygon_rejects_self_intersection_and_accepts_closed_vertex_list():
    with pytest.raises(ValueError, match="self-intersect"):
        Polygon(((0, 0), (2, 2), (0, 2), (2, 0)))
    square = Polygon(((0, 0), (1, 0), (1, 1), (0, 1), (0, 0)))
    assert len(square.vertices) == 4


def test_boolean_region_validates_arity_and_ambient_dimension():
    from semialg.standard_regions import BooleanRegion

    with pytest.raises(ValueError, match="exactly two"):
        BooleanRegion("difference", (Interval(0, 1),))
    with pytest.raises(ValueError, match="exactly one"):
        BooleanRegion("complement", (Interval(0, 1), Interval(2, 3)))
    with pytest.raises(ValueError, match="ambient dimension"):
        BooleanRegion("union", (Interval(0, 1), Ball((0, 0), 1)))


def test_boolean_intersection_dimension_handles_dimension_drop():
    assert BooleanRegion.intersection(Interval(0, 1), Interval(1, 2)).dimension() == 0


def test_concave_polygon_integration_uses_nonoverlapping_triangulation():
    x, y = sp.symbols("x y", real=True)
    polygon = Polygon(((0, 0), (3, 0), (3, 3), (2, 3), (2, 1), (1, 1), (1, 3), (0, 3)))
    assert integrate_over_standard_region(1, polygon, (x, y)) == 7


def test_empty_and_open_degenerate_regions_have_empty_set_dimension():
    from semialg import FinitePointSet

    assert FinitePointSet(()).dimension() == -1
    assert Interval(0, 0).dimension() == 0
    assert Interval(0, 0, lower_closed=False).dimension() == -1
    with pytest.raises(ValueError, match="at least one vertex"):
        Simplex(())
    degenerate = Simplex(((0, 0), (1, 0), (0, 1), (1, 1)))
    assert degenerate.dimension() == 2
    assert not degenerate.is_nondegenerate()


def test_standard_region_dimension_never_exceeds_ambient_dimension():
    regions = (
        Box(((0, 0), (0, 1))),
        Simplex(((0, 0), (1, 0), (2, 0))),
        Parallelogram((0, 0), ((1, 0), (2, 0))),
        Ball((0, 0), 0),
        Sphere((0, 0, 0), 0),
        Cylinder((0, 0, 0), (0, 0, 1), 0),
        Stadium((0, 0), (1, 0), 0),
        Capsule((0, 0, 0), (0, 0, 1), 0),
    )
    assert all(region.dimension() <= region.ambient_dimension() for region in regions)


def test_transformed_region_validates_base_and_coordinate_count_early():
    from semialg import TransformedRegion

    base = Box(((0, 1), (0, 1)))
    with pytest.raises(ValueError, match="base variable count"):
        TransformedRegion(base, (0,), (sp.Symbol("u"),))
    with pytest.raises(TypeError, match="StandardRegion"):
        TransformedRegion(sp.true, (0,), ())  # type: ignore[arg-type]


def test_boolean_region_rejects_nonregion_members_early():
    from semialg.standard_regions import BooleanRegion

    with pytest.raises(TypeError, match="StandardRegion"):
        BooleanRegion("union", (Interval(0, 1), sp.true))  # type: ignore[arg-type]
