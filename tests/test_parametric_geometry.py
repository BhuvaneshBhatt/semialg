import pytest
import sympy as sp

from semialg.affine_geometry import analyze_affine_map
from semialg.geometry_queries import distance_to_region
from semialg.implicit_geometry import (
    decompose_implicit_formula,
    extract_symbolic_box_bounds,
    semialgebraic_level_function,
)
from semialg.map_degree import parametric_map_degree
from semialg.parametric_geometry import bounded_parametric_cover
from semialg.polyhedral_clipping import clip_affine_subspace_to_box
from semialg.region_analysis import region_boundary_result
from semialg.regions.operations import region_dimension
from semialg.standard_regions import BoxRegion, ParametricRegion, SimplexRegion


def test_bounded_parametric_cover_reuses_structured_region_geometry():
    box = BoxRegion(((0, 2), (-1, 1)))
    cover = bounded_parametric_cover(box)
    assert cover.exact
    assert len(cover.charts) == 1
    assert cover.certified_dimension() == 2
    assert region_dimension(box) == 2

    triangle = SimplexRegion(((0, 0), (1, 0), (0, 1)))
    tri_cover = bounded_parametric_cover(triangle)
    assert tri_cover.certified_dimension() == 2
    assert region_dimension(triangle) == 2


def test_bounded_parametric_cover_accepts_exact_clipping_box_for_formula():
    x, y = sp.symbols("x y", real=True)
    cover = bounded_parametric_cover(x**2 + y**2 <= 4, (x, y), ((-1, 1), (-2, 2)))
    assert len(cover.charts) == 1
    chart = cover.charts[0]
    assert chart.mapping == (x, y)
    assert chart.bounds == ((x, -1, 1), (y, -2, 2))


def test_parametric_region_chart_certifies_map_rank_dimension():
    t = sp.symbols("t", real=True)
    parabola = ParametricRegion((t,), ((t, -1, 1),), (t, t**2))
    assert bounded_parametric_cover(parabola).certified_dimension() == 1
    assert region_dimension(parabola) == 1


def test_boundary_result_distinguishes_included_and_excluded_boundary_cells():
    x = sp.symbols("x", real=True)
    result = region_boundary_result(sp.And(x >= 0, x < 1), (x,))
    assert len(result.strata) == 2
    assert len(result.included_strata) == 1
    assert len(result.excluded_strata) == 1
    assert all(stratum.dimension == 0 for stratum in result.strata)
    assert result.cad_result is not None


def test_distance_to_smooth_circle_uses_critical_point_backend():
    x, y = sp.symbols("x y", real=True)
    result = distance_to_region((2, 0), sp.Eq(x**2 + y**2, 1), (x, y), return_result=True)
    assert result.distance == 1
    assert result.optimization.method == "distance_critical_points"
    assert result.optimization.certified
    assert result.nearest_points == ({x: 1, y: 0},)


def test_affine_map_analysis_detects_scaled_isometry_and_inverse():
    x, y = sp.symbols("x y", real=True)
    analysis = analyze_affine_map((2 * x + 3, 2 * y - 5), (x, y))
    assert analysis.rank == 2
    assert analysis.invertible
    assert analysis.scaled_isometry
    assert analysis.scale_factor == 2
    assert not analysis.isometry
    assert analysis.inverse_matrix == sp.ImmutableMatrix(
        [[sp.Rational(1, 2), 0], [0, sp.Rational(1, 2)]]
    )


def test_parametric_map_degree_counts_generic_fibers():
    x, y = sp.symbols("x y", real=True)
    assert parametric_map_degree((x**2,), (x,)).degree == 2
    assert parametric_map_degree((x + y, x - y), (x, y)).degree == 1
    assert parametric_map_degree((x**2, y**2), (x, y)).degree == 4


def test_affine_plane_box_clipping_returns_exact_rectangle_vertices():
    clip = clip_affine_subspace_to_box(
        (0, 0, 0), ((1, 0, 0), (0, 1, 0)), ((-1, 1), (-2, 2), (-3, 3))
    )
    assert clip.dimension == 2
    assert set(clip.vertices) == {(-1, -2, 0), (-1, 2, 0), (1, -2, 0), (1, 2, 0)}


def test_level_function_relaxes_strict_and_handles_reversed_relations():
    x = sp.symbols("x", real=True)
    level = semialgebraic_level_function(sp.And(1 < x, x <= 3), (x,))
    assert sp.simplify(level - sp.Max(1 - x, x - 3)) == 0
    assert semialgebraic_level_function(sp.Or(x > 2, x < -1), (x,)) == sp.Min(2 - x, x + 1)
    with pytest.raises(NotImplementedError):
        semialgebraic_level_function(sp.Eq(x, 0), (x,))


def test_implicit_decomposition_handles_disjunction_strictness_and_reversed_orientation():
    x = sp.symbols("x", real=True)
    pieces = decompose_implicit_formula(
        sp.Or(sp.And(x > 1, x != 2), sp.And(3 >= x, sp.Eq(x, 0))), (x,)
    )
    assert len(pieces) == 2
    assert any(piece.inequalities == (1 - x,) and piece.equalities == () for piece in pieces)
    assert any(piece.inequalities == (x - 3,) and piece.equalities == (x,) for piece in pieces)


def test_symbolic_box_bounds_merge_redundant_and_symbolic_reversed_bounds():
    x, y, a, b = sp.symbols("x y a b", real=True)
    bounds = extract_symbolic_box_bounds(
        sp.And(a <= x, x >= a - 1, x <= b, b >= x, -2 < y, y < 5),
        (x, y),
    )
    assert bounds.limits[0] == (x, a, b)
    assert bounds.limits[1] == (y, -2, 5)
    with pytest.raises(NotImplementedError):
        extract_symbolic_box_bounds(sp.Or(x < 0, x > 1), (x,))


def test_affine_analysis_does_not_overclaim_symbolic_invertibility():
    x, a = sp.symbols("x a", real=True)
    analysis = analyze_affine_map((a * x,), (x,))
    assert analysis.rank == 1  # generic symbolic rank
    assert analysis.invertible is None
    assert analysis.inverse_matrix is None


def test_distance_critical_backend_handles_smooth_complete_intersection():
    x, y, z = sp.symbols("x y z", real=True)
    variety = sp.And(sp.Eq(x, 0), sp.Eq(y, 0))
    result = distance_to_region((2, 3, 5), variety, (x, y, z), return_result=True)
    assert result.squared_distance == 13
    assert result.optimization.method == "distance_critical_points"
    assert result.nearest_points == ({x: 0, y: 0, z: 5},)


def test_parametric_region_exposes_generic_map_degree():
    t = sp.symbols("t", real=True)
    region = ParametricRegion((t,), ((t, -1, 1),), (t**2,))
    assert region.generic_map_degree().degree == 2
