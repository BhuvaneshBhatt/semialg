import pytest
import sympy as sp

from semialg import Ball, Ellipsoid, Sphere
from semialg.standard_regions import BallRegion, SphereRegion


def test_ball_is_arbitrary_dimensional_closed_body():
    x, y, z = sp.symbols("x y z", real=True)
    ball = Ball((1, -1, 2), 3)
    assert ball.dimension() == 3
    assert ball.ambient_dimension() == 3
    assert ball.affine_hull().dimension() == 3
    assert ball.contains((1, -1, 5))
    assert not ball.contains((1, -1, 6))
    assert (
        sp.simplify(ball.as_formula((x, y, z)).lhs - ((x - 1) ** 2 + (y + 1) ** 2 + (z - 2) ** 2))
        == 0
    )
    assert ball.shape_matrix == sp.eye(3) * 9


def test_sphere_is_boundary_not_filled_ball():
    sphere = Sphere((0, 0, 0), 2)
    assert sphere.dimension() == 2
    assert sphere.ambient_dimension() == 3
    assert sphere.contains((0, 0, 2))
    assert not sphere.contains((0, 0, 0))
    assert sphere.affine_hull().dimension() == 3


def test_zero_radius_ball_and_sphere_have_point_affine_hull():
    for region in (Ball((2, 3), 0), Sphere((2, 3), 0)):
        assert region.dimension() == 0
        hull = region.affine_hull()
        assert hull.dimension() == 0
        assert hull.point == (2, 3)


def test_ellipsoid_matrix_formula_and_affine_hull():
    x, y = sp.symbols("x y", real=True)
    ellipsoid = Ellipsoid((1, 2), sp.diag(4, 9))
    assert ellipsoid.dimension() == 2
    assert ellipsoid.affine_hull().dimension() == 2
    assert ellipsoid.contains((3, 2))
    assert not ellipsoid.contains((4, 2))
    formula = ellipsoid.as_formula((x, y))
    expected = (x - 1) ** 2 / 4 + (y - 2) ** 2 / 9 <= 1
    assert sp.simplify(formula.lhs - expected.lhs) == 0


def test_ellipsoid_from_radii_is_canonical_matrix_form():
    ellipsoid = Ellipsoid.from_radii((0, 0, 0), (2, 3, 4))
    assert ellipsoid.shape_matrix == sp.diag(4, 9, 16)
    assert ellipsoid.construction_conditions == ()


def test_ellipsoid_validation_rejects_invalid_quadrics():
    with pytest.raises(ValueError, match="2 x 2"):
        Ellipsoid((0, 0), sp.eye(3))
    with pytest.raises(ValueError, match="symmetric"):
        Ellipsoid((0, 0), [[1, 1], [0, 1]])
    with pytest.raises(ValueError, match="positive definite"):
        Ellipsoid((0, 0), sp.diag(1, -1))
    with pytest.raises(ValueError, match="positive"):
        Ellipsoid.from_radii((0, 0), (1, 0))


def test_symbolic_ellipsoid_keeps_positive_definiteness_conditions():
    a, b = sp.symbols("a b", real=True)
    ellipsoid = Ellipsoid((0, 0), sp.diag(a, b))
    assert ellipsoid.construction_conditions == (a > 0, a * b > 0)


def test_canonical_bridges_from_existing_radial_regions():
    old_ball = BallRegion((1, 2), 3)
    old_sphere = SphereRegion((1, 2, 3), 4)
    assert Ball.from_region(old_ball) == Ball((1, 2), 3)
    assert Sphere.from_region(old_sphere) == Sphere((1, 2, 3), 4)
    with pytest.raises(TypeError, match="sphere boundary"):
        Ball.from_region(SphereRegion((0, 0), 1))
