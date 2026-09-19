import sympy as sp

from semialg import (
    Ball,
    Hyperplane,
    Line,
    Polytope,
    Ray,
    is_disjoint,
    is_subset,
    minkowski_sum,
    region_intersection,
)
from semialg.standard_regions import TransformedRegion


def square(a=0, b=2):
    return Polytope(((a, a), (b, a), (b, b), (a, b)))


def test_line_clips_polytope_without_cad():
    result = region_intersection(Line((-1, 1), (1, 0)), square())
    assert isinstance(result, TransformedRegion)
    assert result.contains((0, 1))
    assert result.contains((2, 1))
    assert not result.contains((3, 1))


def test_ray_clips_polytope():
    result = region_intersection(Ray((1, 1), (1, 0)), square())
    assert isinstance(result, TransformedRegion)
    assert result.contains((1, 1))
    assert result.contains((2, 1))
    assert not result.contains((0, 1))


def test_line_clips_ball_quadratically():
    result = region_intersection(Line((-3, 0), (1, 0)), Ball((0, 0), 2))
    assert isinstance(result, TransformedRegion)
    assert result.contains((-2, 0))
    assert result.contains((2, 0))
    assert not result.contains((sp.Rational(5, 2), 0))


def test_hyperplane_slice_of_cube_is_polytope():
    cube = Polytope(
        ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
    )
    result = region_intersection(Hyperplane((0, 0, 1), (0, 0, sp.Rational(1, 2))), cube)
    assert isinstance(result, Polytope)
    assert result.dimension() == 2
    assert set(result.vertices) == {
        (0, 0, sp.Rational(1, 2)),
        (1, 0, sp.Rational(1, 2)),
        (1, 1, sp.Rational(1, 2)),
        (0, 1, sp.Rational(1, 2)),
    }


def test_polytope_subset_uses_vertex_convexity():
    inner = Polytope(((0, 0), (1, 0), (1, 1), (0, 1)))
    assert is_subset(inner, square())
    assert not is_subset(square(), inner)


def test_polytope_disjoint_uses_linear_feasibility():
    assert is_disjoint(square(), Polytope(((3, 3), (4, 3), (4, 4), (3, 4))))


def test_polytope_minkowski_sum_from_vertex_sums():
    left = Polytope(((0, 0), (1, 0), (1, 1), (0, 1)))
    right = Polytope(((0, 0), (2, 0), (2, 1), (0, 1)))
    result = minkowski_sum(left, right)
    assert isinstance(result, Polytope)
    assert {(0, 0), (3, 0), (3, 2), (0, 2)}.issubset(set(result.vertices))
    assert result.dimension() == 2


def test_hyperplane_slice_of_ball_is_intrinsic_ball():
    result = region_intersection(Hyperplane((0, 0, 1), (0, 0, 1)), Ball((0, 0, 0), 2))
    assert isinstance(result, TransformedRegion)
    assert result.dimension() == 2
    assert result.contains((0, 0, 1))
    assert result.contains((sp.sqrt(3), 0, 1))
    assert not result.contains((2, 0, 1))
