import pytest
import sympy as sp

from semialg import (
    Ball,
    BooleanRegion,
    Capsule,
    Cone,
    Cylinder,
    Interval,
    Parallelogram,
    ParametricRegion,
    Polygon,
    Simplex,
    Sphere,
    Stadium,
    integrate_over_region,
)
from semialg.parametric_integration import (
    integrate_over_parametric_region,
    metric_jacobian_factor,
)

pytestmark = pytest.mark.slow


def test_polygon_tetrahedron_and_parallelogram_integrals():
    x, y, z = sp.symbols("x y z", real=True)
    square = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    assert sp.simplify(integrate_over_region(1, square, [x, y]) - 1) == 0
    assert sp.simplify(integrate_over_region(x + y, square, [x, y]) - 1) == 0

    tet = Simplex([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
    assert sp.simplify(integrate_over_region(1, tet, [x, y, z]) - sp.Rational(1, 6)) == 0

    para = Parallelogram((0, 0), [(2, 0), (0, 3)])
    assert sp.simplify(integrate_over_region(1, para, [x, y]) - 6) == 0


def test_ball_sphere_shell_cylinder_cone_capsule_constant_measures():
    x, y, z = sp.symbols("x y z", real=True)
    assert (
        sp.simplify(
            integrate_over_region(1, Ball((0, 0, 0), 2), [x, y, z]) - sp.Rational(32, 3) * sp.pi
        )
        == 0
    )
    assert sp.simplify(integrate_over_region(1, Sphere((0, 0, 0), 2), [x, y, z]) - 16 * sp.pi) == 0
    assert (
        sp.simplify(
            integrate_over_region(1, Cylinder((0, 0, 0), (0, 0, 3), 2), [x, y, z]) - 12 * sp.pi
        )
        == 0
    )
    assert (
        sp.simplify(integrate_over_region(1, Cone((0, 0, 0), (0, 0, 3), 2), [x, y, z]) - 4 * sp.pi)
        == 0
    )
    assert (
        sp.simplify(integrate_over_region(1, Stadium((0, 0), (2, 0), 1), [x, y]) - (4 + sp.pi)) == 0
    )
    assert (
        sp.simplify(
            integrate_over_region(1, Capsule((0, 0, 0), (0, 0, 3), 1), [x, y, z])
            - (3 * sp.pi + sp.Rational(4, 3) * sp.pi)
        )
        == 0
    )


def test_boolean_symmetric_difference_and_parametric_metric_jacobian():
    x = sp.symbols("x", real=True)
    a = Interval(0, 2)
    b = Interval(1, 3)
    assert (
        sp.simplify(integrate_over_region(1, BooleanRegion.symmetric_difference(a, b), [x]) - 2)
        == 0
    )
    assert sp.simplify(integrate_over_region(1, BooleanRegion.union(a, b), [x]) - 3) == 0

    t, _u, X, Y = sp.symbols("t u X Y", real=True)
    segment = ParametricRegion([t], [(t, 0, 1)], [t, t])
    assert (
        sp.simplify(metric_jacobian_factor(segment.mapping, segment.parameters) - sp.sqrt(2)) == 0
    )
    assert sp.simplify(integrate_over_parametric_region(1, [X, Y], segment) - sp.sqrt(2)) == 0
