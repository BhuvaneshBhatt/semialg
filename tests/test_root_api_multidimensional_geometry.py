"""Metamorphic/boundary dimensions for multidimensional geometry APIs."""

from __future__ import annotations

import sympy as sp
from hypothesis import given, settings
from hypothesis import strategies as st

import semialg as sa

X, Y, Z = sp.symbols("x y z", real=True)


def test_dodecahedron_is_translation_covariant_and_has_twenty_vertices():
    a = sa.Dodecahedron(center=(0, 0, 0), edge=2)
    b = sa.Dodecahedron(center=(3, -2, 5), edge=2)
    assert len(a.vertices) == len(b.vertices) == 20
    assert {
        tuple(sp.simplify(q - p) for p, q in zip(u, v, strict=True))
        for u, v in zip(a.vertices, b.vertices, strict=True)
    } == {(3, -2, 5)}


def test_active_constraints_change_exactly_at_triangle_boundary():
    region = (X >= 0) & (Y >= 0) & (X + Y <= 1)
    interior = sa.active_constraints(region, (sp.Rational(1, 3), sp.Rational(1, 3)), (X, Y))
    vertex = sa.active_constraints(region, (0, 0), (X, Y))
    assert interior.constraints == ()
    assert len(vertex.constraints) == 2


def test_bounding_box_is_translation_covariant():
    region = (X >= -2) & (X <= 3) & (Y >= 1) & (Y <= 4)
    box = sa.bounding_box(region, (X, Y))
    moved = sa.bounding_box(sa.translate(region, (7, -5), (X, Y)), (X, Y))
    assert box == {X: (-2, 3), Y: (1, 4)}
    assert moved == {X: (5, 10), Y: (-4, -1)}


def test_canonicalize_polygon_is_invariant_under_cyclic_vertex_rotation():
    p = sa.Polygon(((0, 0), (2, 0), (2, 1), (0, 1)))
    q = sa.Polygon(((2, 1), (0, 1), (0, 0), (2, 0)))
    assert sa.canonicalize_polygon(p) == sa.canonicalize_polygon(q)


def test_canonicalize_region_is_idempotent_on_polygon():
    region = sa.Polygon(((0, 0), (2, 0), (0, 1)))
    canonical = sa.canonicalize_region(region)
    assert sa.canonicalize_region(canonical) == canonical


def test_closest_points_and_distance_between_regions_share_the_same_optimum():
    left = X <= 0
    right = X >= 3
    pairs = sa.closest_points(left, right, (X,))
    assert sa.distance_between_regions(left, right, (X,)) == 3
    assert pairs and pairs[0][0][X] == 0 and pairs[0][1][X] == 3


def test_regularization_operators_are_idempotent_on_closed_and_open_intervals():
    closed = (X >= 0) & (X <= 1)
    opened = (X > 0) & (X < 1)
    assert sa.is_equal(sa.closure_of_interior(closed, (X,)), closed, (X,))
    assert sa.is_equal(sa.interior_of_closure(opened, (X,)), opened, (X,))
    assert sa.is_regular_closed_region(closed, (X,))
    assert sa.is_regular_open_region(opened, (X,))


def test_component_constraint_descriptions_cover_each_axis_component():
    descriptions = sa.component_constraint_descriptions(sp.Eq(X * Y, 0), (X, Y))
    assert len(descriptions) == 2
    assert any(
        sa.contains_point(item.component_formula, {X: 0, Y: 2}, (X, Y)) for item in descriptions
    )
    assert any(
        sa.contains_point(item.component_formula, {X: 2, Y: 0}, (X, Y)) for item in descriptions
    )


def test_critical_value_image_contains_stationary_value_on_compact_interval():
    image = sa.critical_value_image(X**2, (X >= -2) & (X <= 3), (X,))
    assert image is not None
    assert 0 in image.isolated_values


def test_polytope_decomposition_preserves_triangle_area():
    poly = sa.Polytope(((0, 0), (2, 0), (0, 2)))
    result = sa.decompose_polytope(poly)
    assert result.certified and result.conforming
    assert result.simplices == ((0, 1, 2),)


def test_deduplicate_indexed_vertices_is_idempotent_and_reindexes_cells():
    first = sa.deduplicate_indexed_vertices(((0, 0), (1, 0), (0, 0)), ((0, 1), (2, 1)))
    second = sa.deduplicate_indexed_vertices(first.vertices, first.cells)
    assert (first.vertices, first.cells) == (second.vertices, second.cells)
    assert first.cells == ((0, 1), (0, 1))


def test_distance_set_and_squared_distance_range_agree_under_squaring():
    distances = sa.distance_set((X >= 0) & (X <= 2), variables=(X,))
    squared = sa.squared_distance_range((X >= 0) & (X <= 2), variables=(X,))
    d = next(s for s in distances.free_symbols if s.name == "d")
    t = next(s for s in squared.free_symbols if s.name == "t")
    for value in (0, 1, 2):
        assert bool(distances.subs(d, value)) == bool(squared.subs(t, value**2))


def test_euler_characteristic_is_additive_for_two_disjoint_closed_intervals():
    region = ((X >= -2) & (X <= -1)) | ((X >= 1) & (X <= 3))
    assert sa.euler_characteristic(region, (X,)) == 2


def test_extrema_set_is_union_of_argmin_and_argmax_sets():
    region = (X >= -2) & (X <= 3)
    extrema = sa.extrema_set(X**2, region, (X,))
    expected = sa.region_union(sa.argmin_set(X**2, region, (X,)), sa.argmax_set(X**2, region, (X,)))
    assert sa.is_equal(extrema, expected, (X,))


def test_fiber_substitution_matches_direct_substitution():
    region = (X**2 + Y**2 <= 4) & (X >= 0)
    fib = sa.fiber(region, {X: 1})
    assert sp.simplify_logic(fib ^ (Y**2 <= 3)) is sp.false


def test_geodesic_refinement_preserves_sphere_radius():
    mesh = sa.geodesic_refinement(sa.Octahedron(), levels=0, center=(0, 0, 0), radius=2)
    assert all(sp.simplify(sum(c * c for c in vertex) - 4) == 0 for vertex in mesh.vertices)


def test_has_empty_interior_distinguishes_line_from_strip():
    assert sa.has_empty_interior(sp.Eq(Y, 0), (X, Y)) is True
    assert sa.has_empty_interior((Y >= 0) & (Y <= 1), (X, Y)) is False


def test_implied_and_redundant_polynomial_inequality_agree():
    region = (X >= 0) & (2 * X >= 0) & (X <= 1)
    implication = sa.implied_polynomial_inequality(region, 2 * X >= 0, (X,))
    redundant = sa.redundant_polynomial_inequalities(region, (X,))
    assert implication.certified and redundant


def test_moment_and_inertia_tensor_obey_parallel_axis_symmetry_at_origin():
    square = (X >= -1) & (X <= 1) & (Y >= -1) & (Y <= 1)
    assert sa.moment_matrix(square, (X, Y)) == sp.diag(sp.Rational(1, 3), sp.Rational(1, 3))
    assert sa.inertia_tensor(square, (X, Y)) == sp.diag(sp.Rational(4, 3), sp.Rational(4, 3))


def _polygonal_with_hole():
    component = sa.PolygonalComponent(
        ((0, 0), (0, 4), (4, 4), (4, 0)), (((1, 1), (2, 1), (2, 2), (1, 2)),)
    )
    return sa.PolygonalSet((component,))


def test_polygon_shell_accessors_partition_outer_and_inner_boundaries():
    region = _polygonal_with_hole()
    assert len(sa.outer_polygons(region)) == 1
    assert len(sa.inner_polygons(region)) == 1
    assert len(sa.polygon_vertices(region)) == 8


def _cube_boundary():
    cube = sa.Hexahedron(
        ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
    )
    return sa.convert_region(cube, "boundary")


def test_polyhedral_shell_accessors_preserve_outer_shell_without_cavities():
    boundary = _cube_boundary()
    assert len(sa.outer_polyhedra(boundary)) == 1
    assert sa.inner_polyhedra(boundary) == ()
    assert len(sa.polyhedron_vertices(boundary)) == 8
    assert len(sa.polyhedron_face_indices(boundary)) == 1


def test_topological_predicates_cover_open_closed_dense_and_full_dimensional_cases():
    closed = (X >= 0) & (X <= 1)
    opened = (X > 0) & (X < 1)
    assert sa.is_closed(closed, (X,)) and not sa.is_open(closed, (X,))
    assert sa.is_open(opened, (X,)) and not sa.is_closed(opened, (X,))
    assert sa.is_dense_in(sp.Ne(X, 0), sp.true, (X,))
    assert sa.is_full_dimensional(closed, (X,))
    assert sa.is_path_connected(closed, (X,))
    assert sa.intersects(closed, X >= sp.Rational(1, 2), (X,))


def test_smoothness_and_singularity_predicates_agree_on_circle_and_cusp():
    circle = X**2 + Y**2 - 1
    cusp = Y**2 - X**3
    assert sa.is_smooth((circle,), (X, Y))
    assert not sa.is_singular((circle,), {X: 1, Y: 0}, (X, Y))
    assert sa.is_singular((cusp,), {X: 0, Y: 0}, (X, Y))


def test_level_sublevel_superlevel_are_nested_at_same_value():
    level = sa.level_set(X**2, 1, sp.true)
    sub = sa.sublevel_set(X**2, 1, sp.true)
    sup = sa.superlevel_set(X**2, 1, sp.true)
    assert sa.is_subset(level, sub, (X,)) and sa.is_subset(level, sup, (X,))


def test_linear_image_commutes_with_positive_scalar_interval_bounds():
    interval = (X >= -1) & (X <= 2)
    image = sa.linear_image(interval, [[3]], (X,))
    assert sa.is_equal(image, (X >= -3) & (X <= 6), (X,))


def test_nearest_point_realizes_distance_to_region():
    region = (X >= 0) & (X <= 2)
    nearest = sa.nearest_point((5,), region, (X,))
    assert nearest and nearest[0][X] == 2
    assert sa.distance_to_region((5,), region, (X,)) == 3


def test_parameterization_critical_locus_and_values_for_parabola():
    t = sp.Symbol("t", real=True)
    u, v = sp.symbols("u v", real=True)
    locus = sa.parameterization_critical_locus((t**2, t**3), sp.true, (t,))
    values = sa.parameterization_critical_values(
        (t**2, t**3), sp.true, (t,), image_variables=(u, v)
    )
    assert bool(locus.subs(t, 0))
    assert bool(values.subs({u: 0, v: 0}))


def test_polyhedral_intersection_and_boolean_intersection_agree():
    a = sa.Polytope(((0, 0), (2, 0), (2, 2), (0, 2)))
    b = sa.Polytope(((1, 1), (3, 1), (3, 3), (1, 3)))
    direct = sa.polyhedral_intersection(a, b)
    boolean = sa.polyhedral_boolean("intersection", a, b)
    assert sa.canonicalize_region(direct) == sa.canonicalize_region(boolean)


def test_polynomial_constraints_are_invariant_under_conjunct_reordering():
    a = sa.polynomial_constraints((X >= 0) & (Y >= 0) & (X + Y <= 1), (X, Y))
    b = sa.polynomial_constraints((X + Y <= 1) & (Y >= 0) & (X >= 0), (X, Y))
    assert a == b


@given(st.integers(3, 8), st.integers(1, 20))
@settings(max_examples=12, deadline=None)
def test_random_polygon_is_seed_reproducible_and_canonicalizable(vertex_count, seed):
    a = sa.random_polygon(vertex_count=vertex_count, seed=seed)
    b = sa.random_polygon(vertex_count=vertex_count, seed=seed)
    assert a == b
    assert sa.canonicalize_polygon(a) == sa.canonicalize_polygon(b)


def test_random_polytope_is_seed_reproducible_and_full_dimensional():
    a = sa.random_polytope(dimension=2, point_count=6, seed=7)
    b = sa.random_polytope(dimension=2, point_count=6, seed=7)
    assert a == b
    assert a.dimension() == 2


def test_boundary_and_singular_result_objects_agree_with_regular_locus():
    disk = X**2 + Y**2 <= 1
    boundary = sa.region_boundary_result(disk, (X, Y))
    singular = sa.region_singular_locus_result(disk, (X, Y))
    regular = sa.region_regular_locus(disk, (X, Y))
    assert boundary.formula is not None
    assert singular.formula is not None
    assert regular is not None


def test_symmetric_difference_is_commutative():
    a = X <= 0
    b = X >= 1
    assert sa.is_equal(
        sa.region_symmetric_difference(a, b), sa.region_symmetric_difference(b, a), (X,)
    )


def test_region_variables_ignore_parameters_when_coordinates_are_explicit():
    a = sp.Symbol("a", real=True)
    assert sa.region_variables(X >= a, (X,), kind="coordinates") == (X,)


def test_relative_boundary_and_interior_partition_closed_segment_in_affine_hull():
    segment = sp.Eq(Y, 0) & (X >= 0) & (X <= 1)
    interior = sa.relative_interior(segment, (X, Y))
    boundary = sa.relative_boundary(segment, (X, Y))
    assert sa.contains_point(interior, {X: sp.Rational(1, 2), Y: 0}, (X, Y))
    assert sa.contains_point(boundary, {X: 0, Y: 0}, (X, Y))
    assert not sa.contains_point(interior, {X: 0, Y: 0}, (X, Y))


def test_semialgebraic_tangent_cone_of_halfline_is_itself():
    cone = sa.semialgebraic_tangent_cone(X >= 0, (0,), (X,))
    d = next(iter(cone.free_symbols))
    assert bool(cone.subs(d, 1)) and not bool(cone.subs(d, -1))


def test_subdivide_faces_growth_and_tetrahedralize_cells_preserve_index_domain():
    refined = sa.subdivide_triangular_faces(sa.Octahedron(), levels=1)
    assert len(refined.faces) == 32
    tetra = sa.tetrahedralize_cells(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)), ((0, 1, 2, 3),))
    assert tetra.tetrahedra == ((0, 1, 2, 3),)


def test_tangent_dimension_matches_line_and_plane_dimensions():
    assert sa.tangent_dimension((Y,), {X: 0, Y: 0}, (X, Y)) == 1
    assert sa.tangent_dimension((Z,), {X: 0, Y: 0, Z: 0}, (X, Y, Z)) == 2
