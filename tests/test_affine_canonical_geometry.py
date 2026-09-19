import pytest
import sympy as sp

from semialg import AffineHalfSpace, AffineSpace, HalfSpace, Hyperplane, Line, Point, Ray


def test_point_geometry_and_formula():
    x, y = sp.symbols("x y", real=True)
    point = Point((1, 2))
    assert point.coordinates == (sp.Integer(1), sp.Integer(2))
    assert point.dimension() == 0
    assert point.ambient_dimension() == 2
    assert point.as_formula((x, y)) == sp.And(sp.Eq(x, 1), sp.Eq(y, 2))
    assert point.affine_hull() == AffineSpace((1, 2), ())


def test_affine_space_uses_independent_directions():
    x, y, z = sp.symbols("x y z", real=True)
    space = AffineSpace((1, 2, 3), ((1, 0, 0), (2, 0, 0), (0, 1, 0)))
    assert space.dimension() == 2
    assert space.ambient_dimension() == 3
    assert len(space.directions) == 2
    assert space.as_formula((x, y, z), eliminate=True) == sp.Eq(z, 3)


def test_hyperplane_and_halfspace_lower_directly():
    x, y = sp.symbols("x y", real=True)
    plane = Hyperplane((2, -1), (3, 4))
    half = HalfSpace((2, -1), (3, 4))
    linear = 2 * (x - 3) - (y - 4)
    assert plane.as_formula((x, y)) == sp.Eq(linear, 0)
    assert half.as_formula((x, y)) == (linear <= 0)
    assert plane.dimension() == 1
    assert half.dimension() == 2
    assert plane.affine_hull().dimension() == 1
    assert half.affine_hull().dimension() == 2


def test_line_and_ray_formulas():
    x, y = sp.symbols("x y", real=True)
    line = Line((1, 2), (2, 3))
    ray = Ray((1, 2), (2, 3))
    assert line.dimension() == ray.dimension() == 1
    assert line.affine_hull() == line
    assert ray.affine_hull() == line
    assert line.contains((5, 8))
    assert line.contains((-1, -1))
    assert ray.contains((5, 8))
    assert not ray.contains((-1, -1))
    assert line.as_formula((x, y), eliminate=True).subs({x: 5, y: 8}) is sp.true


def test_affine_halfspace_intrinsic_geometry():
    x, y, z = sp.symbols("x y z", real=True)
    region = AffineHalfSpace((0, 0, 1), ((1, 0, 0),), (0, 1, 0))
    assert region.dimension() == 2
    assert region.ambient_dimension() == 3
    assert region.affine_hull() == AffineSpace((0, 0, 1), ((1, 0, 0), (0, 1, 0)))
    assert region.contains((4, 3, 1))
    assert not region.contains((4, -3, 1))
    formula = region.as_formula((x, y, z), eliminate=True)
    for values, expected in [
        ({x: 4, y: 3, z: 1}, True),
        ({x: 4, y: -3, z: 1}, False),
        ({x: 0, y: 0, z: 1}, True),
        ({x: 0, y: 2, z: 0}, False),
    ]:
        assert bool(formula.subs(values)) is expected


def test_invalid_affine_geometry_is_rejected():
    with pytest.raises(ValueError, match="nonzero"):
        Hyperplane((0, 0), (0, 0))
    with pytest.raises(ValueError, match="nonzero"):
        Line((0, 0), (0, 0))
    with pytest.raises(ValueError, match="boundary span"):
        AffineHalfSpace((0, 0), ((1, 0),), (2, 0))
    with pytest.raises(ValueError, match="same dimension"):
        HalfSpace((1, 0), (0, 0, 0))
