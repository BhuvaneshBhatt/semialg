import pytest
import sympy as sp

from semialg import (
    BallRegion,
    CapsuleRegion,
    ConeRegion,
    CylinderRegion,
    IntervalRegion,
    ParallelogramRegion,
    ParametricRegion,
    PolygonRegion,
    RegionSymmetricDifference,
    RegionUnion,
    SphereRegion,
    StadiumRegion,
    TetrahedronRegion,
    integrate_over_parametric_region,
    integrate_over_region,
    metric_jacobian_factor,
)

pytestmark = pytest.mark.slow


def test_polygon_tetrahedron_and_parallelogram_integrals():
    x, y, z = sp.symbols("x y z", real=True)
    square = PolygonRegion([(0, 0), (1, 0), (1, 1), (0, 1)])
    assert sp.simplify(integrate_over_region(1, square, [x, y]) - 1) == 0
    assert sp.simplify(integrate_over_region(x + y, square, [x, y]) - 1) == 0

    tet = TetrahedronRegion([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
    assert sp.simplify(integrate_over_region(1, tet, [x, y, z]) - sp.Rational(1, 6)) == 0

    para = ParallelogramRegion((0, 0), [(2, 0), (0, 3)])
    assert sp.simplify(integrate_over_region(1, para, [x, y]) - 6) == 0


def test_ball_sphere_shell_cylinder_cone_capsule_constant_measures():
    x, y, z = sp.symbols("x y z", real=True)
    assert (
        sp.simplify(
            integrate_over_region(1, BallRegion((0, 0, 0), 2), [x, y, z])
            - sp.Rational(32, 3) * sp.pi
        )
        == 0
    )
    assert (
        sp.simplify(integrate_over_region(1, SphereRegion((0, 0, 0), 2), [x, y, z]) - 16 * sp.pi)
        == 0
    )
    assert (
        sp.simplify(
            integrate_over_region(1, CylinderRegion((0, 0, 0), (0, 0, 3), 2), [x, y, z])
            - 12 * sp.pi
        )
        == 0
    )
    assert (
        sp.simplify(
            integrate_over_region(1, ConeRegion((0, 0, 0), (0, 0, 3), 2), [x, y, z]) - 4 * sp.pi
        )
        == 0
    )
    assert (
        sp.simplify(
            integrate_over_region(1, StadiumRegion((0, 0), (2, 0), 1), [x, y]) - (4 + sp.pi)
        )
        == 0
    )
    assert (
        sp.simplify(
            integrate_over_region(1, CapsuleRegion((0, 0, 0), (0, 0, 3), 1), [x, y, z])
            - (3 * sp.pi + sp.Rational(4, 3) * sp.pi)
        )
        == 0
    )


def test_boolean_symmetric_difference_and_parametric_metric_jacobian():
    x = sp.symbols("x", real=True)
    a = IntervalRegion(0, 2)
    b = IntervalRegion(1, 3)
    assert sp.simplify(integrate_over_region(1, RegionSymmetricDifference(a, b), [x]) - 2) == 0
    assert sp.simplify(integrate_over_region(1, RegionUnion(a, b), [x]) - 3) == 0

    t, u, X, Y = sp.symbols("t u X Y", real=True)
    segment = ParametricRegion([t], [(t, 0, 1)], [t, t])
    assert (
        sp.simplify(metric_jacobian_factor(segment.mapping, segment.parameters) - sp.sqrt(2)) == 0
    )
    assert sp.simplify(integrate_over_parametric_region(1, [X, Y], segment) - sp.sqrt(2)) == 0
