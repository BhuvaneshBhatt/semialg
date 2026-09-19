import pytest
import sympy as sp

from semialg import Polygon, RegularPolygon, Simplex, Triangle


def _distance(a, b):
    return sp.sqrt(sum((x - y) ** 2 for x, y in zip(a, b, strict=True)))


def test_triangle_direct_constructor_returns_simplex():
    region = Triangle(((0, 0), (1, 0), (0, 1)))
    assert type(region) is Simplex
    assert region.dimension() == 2


def test_triangle_from_sides_matches_all_three_lengths():
    region = Triangle.from_sides(3, 4, 5)
    a, b, c = region.vertices
    assert sp.simplify(_distance(b, c) - 3) == 0
    assert sp.simplify(_distance(a, c) - 4) == 0
    assert sp.simplify(_distance(a, b) - 5) == 0


def test_triangle_from_sas_matches_sides_and_angle():
    region = Triangle.from_sas(3, sp.pi / 2, 4)
    assert type(region) is Simplex
    assert region.vertices == ((0, 0), (3, 0), (0, 4))


def test_triangle_from_asa_matches_reference_triangle():
    region = Triangle.from_asa(sp.pi / 4, sp.sqrt(2), sp.pi / 4)
    assert sp.simplify(region.vertices[2][0] - sp.sqrt(2) / 2) == 0
    assert sp.simplify(region.vertices[2][1] - sp.sqrt(2) / 2) == 0


def test_triangle_from_aas_obeys_sine_law():
    region = Triangle.from_aas(sp.pi / 6, sp.pi / 3, 2)
    _, b, c = region.vertices
    assert sp.simplify(_distance(b, c) - 2) == 0


def test_triangle_constructors_reject_certified_invalid_data():
    with pytest.raises(ValueError, match="triangle inequalities"):
        Triangle.from_sides(1, 2, 3)
    with pytest.raises(ValueError, match="positive"):
        Triangle.from_sas(0, sp.pi / 3, 1)
    with pytest.raises(ValueError, match="strictly between"):
        Triangle.from_asa(sp.pi / 2, 1, sp.pi / 2)


def test_regular_polygon_returns_canonical_polygon():
    region = RegularPolygon(4)
    assert type(region) is Polygon
    assert region.vertices == ((1, 0), (0, 1), (-1, 0), (0, -1))
    assert region.dimension() == 2


def test_regular_polygon_center_radius_and_rotation():
    region = RegularPolygon(3, center=(2, -1), radius=2, rotation=sp.pi / 2)
    assert region.vertices[0] == (2, 1)
    distances = [sp.simplify(_distance(vertex, (2, -1))) for vertex in region.vertices]
    assert distances == [2, 2, 2]


def test_regular_polygon_rejects_invalid_shape_data():
    with pytest.raises(ValueError, match="at least three"):
        RegularPolygon(2)
    with pytest.raises(ValueError, match="positive"):
        RegularPolygon(4, radius=0)
    with pytest.raises(ValueError, match="two-dimensional"):
        RegularPolygon(4, center=(0, 0, 0))


def test_symbolic_constructors_preserve_validity_conditions():
    a, b, c = sp.symbols("a b c", positive=True)
    r = sp.symbols("r", positive=True)
    triangle = Triangle.from_sides(a, b, c)
    expected = sp.And(a + b > c, a + c > b, b + c > a)
    actual = sp.And(*triangle.construction_conditions)
    assert sp.simplify_logic(sp.Equivalent(actual, expected)) is sp.true
    polygon = RegularPolygon(5, radius=r)
    assert polygon.construction_conditions == ()


def test_formula_lowering_uses_canonical_underlying_types():
    x, y = sp.symbols("x y", real=True)
    triangle = Triangle.from_sides(3, 4, 5)
    square = RegularPolygon(4)
    assert triangle.as_formula((x, y), eliminate=True).subs({x: 1, y: 0}) is not sp.false
    assert square.as_formula((x, y), eliminate=True).subs({x: 0, y: 0}) is not sp.false
