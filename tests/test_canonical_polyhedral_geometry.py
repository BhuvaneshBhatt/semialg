import pytest
import sympy as sp

from semialg import ConicRegion, Parallelepiped, Polygon, Polytope, Simplex
from semialg.standard_regions import ParallelepipedRegion, PolygonRegion, SimplexRegion


def test_polytope_convex_hull_formula_and_affine_hull():
    x, y = sp.symbols("x y", real=True)
    region = Polytope([(0, 0), (2, 0), (0, 1)])
    assert region.dimension() == 2
    assert region.ambient_dimension() == 2
    assert region.contains((sp.Rational(1, 2), sp.Rational(1, 4)))
    assert not region.contains((2, 2))
    assert region.affine_hull().dimension() == 2
    assert region.as_formula((x, y)).has(sp.Symbol)  # structural existential lowering


def test_simplex_requires_affine_independence():
    triangle = Simplex([(0, 0), (1, 0), (0, 1)])
    assert triangle.dimension() == 2
    assert triangle.contains((sp.Rational(1, 4), sp.Rational(1, 4)))
    with pytest.raises(ValueError, match="affinely independent"):
        Simplex([(0, 0), (1, 0), (2, 0)])


def test_polygon_preserves_simple_polygon_semantics():
    polygon = Polygon([(0, 0), (2, 0), (2, 1), (0, 1)])
    assert polygon.dimension() == 2
    assert len(polygon.triangulation()) == 2
    assert polygon.affine_hull().dimension() == 2


def test_conic_region_uses_free_and_nonnegative_generators():
    x, y = sp.symbols("x y", real=True)
    cone = ConicRegion((0, 0), directions=[(1, 0)], rays=[(0, 1)])
    assert cone.dimension() == 2
    assert cone.contains((3, 2))
    assert not cone.contains((3, -2))
    qf = cone.as_formula((x, y), eliminate=True)
    assert sp.simplify_logic(qf.subs({x: 4, y: 3})) is sp.true
    assert sp.simplify_logic(qf.subs({x: 4, y: -3})) is sp.false


def test_parallelepiped_requires_independent_vectors():
    region = Parallelepiped((0, 0), [(2, 0), (0, 3)])
    assert region.dimension() == 2
    assert region.contains((1, 1))
    assert region.affine_hull().dimension() == 2
    with pytest.raises(ValueError, match="linearly independent"):
        Parallelepiped((0, 0), [(1, 0), (2, 0)])


def test_conversion_bridges_preserve_canonical_data():
    assert Simplex.from_region(SimplexRegion([(0, 0), (1, 0), (0, 1)])).vertices == (
        (0, 0),
        (1, 0),
        (0, 1),
    )
    assert Polygon.from_region(PolygonRegion([(0, 0), (1, 0), (1, 1), (0, 1)])).vertices[0] == (
        0,
        0,
    )
    old = ParallelepipedRegion((1, 2), [(1, 0), (0, 2)])
    new = Parallelepiped.from_region(old)
    assert new.origin == old.origin
    assert new.vectors == old.vectors
    assert Polytope.from_region(old).ambient_dimension() == 2
