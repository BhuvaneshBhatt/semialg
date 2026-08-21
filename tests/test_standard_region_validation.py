import pytest
import sympy as sp

from semialg import (
    BallRegion,
    BoxRegion,
    CapsuleRegion,
    CylinderRegion,
    IntervalRegion,
    ParallelepipedRegion,
    ParallelogramRegion,
    ParametricRegion,
    PolygonRegion,
    PrismRegion,
    PyramidRegion,
    RegionDifference,
    RegionIntersection,
    RegionUnion,
    SimplexRegion,
    SphereRegion,
    SphericalShellRegion,
    StadiumRegion,
    integrate_over_parametric_region,
    integrate_over_standard_region,
    reduce_parametric_integral,
)


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
        IntervalRegion(2, 1)
    with pytest.raises(ValueError, match="lower endpoint"):
        BoxRegion(((0, 1), (2, 1)))
    assert integrate_over_standard_region(1, IntervalRegion(1, 2), [x]) == 1


def test_radial_regions_reject_invalid_radii():
    constructors = (
        lambda: BallRegion((0, 0), -1),
        lambda: SphereRegion((0, 0, 0), -1),
        lambda: CylinderRegion((0, 0, 0), (0, 0, 1), -1),
        lambda: StadiumRegion((0, 0), (1, 0), -1),
        lambda: CapsuleRegion((0, 0, 0), (0, 0, 1), -1),
    )
    for make_region in constructors:
        with pytest.raises(ValueError, match="nonnegative"):
            make_region()
    with pytest.raises(ValueError, match="shell radii"):
        SphericalShellRegion((0, 0, 0), (2, 1))
    with pytest.raises(ValueError, match="nonnegative"):
        SphericalShellRegion((0, 0, 0), (-1, 2))


def test_standard_regions_validate_ambient_dimensions():
    with pytest.raises(ValueError, match="same dimension"):
        SimplexRegion(((0, 0), (1, 0, 0)))
    with pytest.raises(ValueError, match="origin dimension"):
        ParallelogramRegion((0, 0), ((1, 0), (0, 1, 0)))
    with pytest.raises(ValueError, match="origin dimension"):
        ParallelepipedRegion((0, 0, 0), ((1, 0, 0), (0, 1)))
    base = PolygonRegion(((0, 0), (1, 0), (0, 1)))
    with pytest.raises(ValueError, match="base ambient dimension"):
        PrismRegion(base, (0, 0, 1))
    with pytest.raises(ValueError, match="base ambient dimension"):
        PyramidRegion(base, (0, 0, 1))
    with pytest.raises(ValueError, match="same dimension"):
        CylinderRegion((0, 0), (0, 0, 1), 1)


def test_interval_measure_identity_for_small_rational_intervals():
    x = sp.Symbol("x", real=True)
    endpoints = tuple(sp.Rational(i, 2) for i in range(-2, 3))
    for lower in endpoints:
        for upper in endpoints:
            if lower <= upper:
                region = IntervalRegion(lower, upper)
                assert integrate_over_standard_region(1, region, [x]) == upper - lower


def test_interval_boolean_measure_identities():
    x = sp.Symbol("x", real=True)
    cases = (
        (IntervalRegion(0, 2), IntervalRegion(1, 3)),
        (IntervalRegion(0, 1), IntervalRegion(2, 3)),
        (IntervalRegion(-1, 3), IntervalRegion(0, 1)),
    )
    for left, right in cases:
        m_left = integrate_over_standard_region(1, left, [x])
        m_right = integrate_over_standard_region(1, right, [x])
        m_inter = integrate_over_standard_region(1, RegionIntersection(left, right), [x])
        m_union = integrate_over_standard_region(1, RegionUnion(left, right), [x])
        assert sp.simplify(m_union + m_inter - m_left - m_right) == 0
        m_diff = integrate_over_standard_region(1, RegionDifference(left, right), [x])
        assert sp.simplify(m_diff - (m_left - m_inter)) == 0


def test_region_difference_of_disjoint_intervals_preserves_left_region():
    x = sp.Symbol("x", real=True)
    left = IntervalRegion(0, 1)
    right = IntervalRegion(2, 3)
    assert integrate_over_standard_region(1, RegionDifference(left, right), [x]) == 1
