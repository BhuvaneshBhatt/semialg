"""Boundary and invariant contracts for thin high-risk root APIs."""

from __future__ import annotations

import sympy as sp

from semialg import (
    Cube,
    Dodecahedron,
    Icosahedron,
    IntervalRegion,
    Octahedron,
    Point,
    Prism,
    Pyramid,
    RegionDifference,
    RegionIntersection,
    RegionSymmetricDifference,
    RegularPolygon,
    SemialgebraicRegion,
    Tetrahedron,
    argmax_set,
    betti_number,
    bounding_box,
    closest_points,
    connected_component_count,
    connected_component_samples,
    contains_point,
    critical_values,
    distance_between_regions,
    distance_set,
    equivalent,
    euler_characteristic,
    extrema_set,
    has_empty_interior,
    is_bounded,
    is_dense_in,
    is_full_dimensional,
    is_path_connected,
    is_regular_closed_region,
    is_regular_open_region,
    is_singular,
    is_smooth,
    level_set,
    linear_image,
    local_dimension,
    nearest_point,
    region_boundary_result,
    region_centroid,
    region_closure_interior,
    region_covariance,
    region_interior_closure,
    region_measure,
    region_product,
    region_regular_locus,
    region_variables,
    semialgebraic_preimage,
    squared_distance_range,
    sublevel_set,
    superlevel_set,
    tangent_dimension,
    topology_summary,
)
from semialg.algebraic import thom_encoding, thom_encodings
from semialg.map_degree import parametric_map_degree
from semialg.parametric_geometry import bounded_parametric_cover
from semialg.roadmaps import roadmap
from semialg.topology.semialgebraic import (
    component_decomposition,
    dimension_strata,
    hardt_trivialization,
    triangulate_region,
)


def test_named_solid_factories_preserve_dimension_and_requested_scale() -> None:
    cube = Cube(center=(1, 2, 3), side=2)
    tetrahedron = Tetrahedron(center=(0, 0, 0), edge=2)
    octahedron = Octahedron(center=(0, 0, 0), edge=2)
    icosahedron = Icosahedron(center=(0, 0, 0), edge=2)
    dodecahedron = Dodecahedron(center=(0, 0, 0), edge=2)
    prism = Prism(((0, 0, 0), (1, 0, 0), (0, 1, 0)), (0, 0, 1))
    pyramid = Pyramid(((0, 0, 0), (1, 0, 0), (0, 1, 0)), (0, 0, 1))

    assert cube.dimension() == 3 and cube.origin == (0, 1, 2)
    assert tetrahedron.dimension() == 3 and len(tetrahedron.vertices) == 4
    assert octahedron.dimension() == 3 and len(octahedron.vertices) == 6
    assert icosahedron.dimension() == 3 and len(icosahedron.vertices) == 12
    assert dodecahedron.dimension() == 3 and len(dodecahedron.vertices) == 20
    assert prism.dimension() == 3 and len(prism.vertices) == 6
    assert pyramid.dimension() == 3 and len(pyramid.vertices) == 4


def test_polygon_and_boolean_region_factories_preserve_membership() -> None:
    left = IntervalRegion(0, 2)
    right = IntervalRegion(1, 3)
    polygon = RegularPolygon(4, radius=1, rotation=sp.pi / 4)
    difference = RegionDifference(left, right)
    intersection = RegionIntersection(left, right)
    symmetric = RegionSymmetricDifference(left, right)

    assert polygon.dimension() == 2 and len(polygon.vertices) == 4
    assert difference.contains((sp.Rational(1, 2),))
    assert not difference.contains((sp.Rational(3, 2),))
    assert intersection.contains((sp.Rational(3, 2),))
    assert not symmetric.contains((sp.Rational(3, 2),))
    assert symmetric.contains((sp.Rational(5, 2),))
    x = sp.Symbol("x", real=True)
    assert contains_point((x >= 0) & (x <= 1), (sp.Rational(1, 2),), (x,))


def test_interval_geometry_queries_cross_check_one_another() -> None:
    x = sp.Symbol("x", real=True)
    interval = (x >= 0) & (x <= 1)
    left_point = sp.Eq(x, 0)
    right_point = sp.Eq(x, 2)

    assert bounding_box(interval, (x,)) == {x: (0, 1)}
    assert closest_points(left_point, right_point, (x,)) == (({x: 0}, {x: 2}),)
    assert distance_between_regions(left_point, right_point, (x,)) == 2
    distance_formula = distance_set(left_point, right_point, (x,))
    distance_symbol = next(iter(distance_formula.free_symbols))
    assert sp.simplify(distance_formula.subs(distance_symbol, 2)) is sp.true
    assert nearest_point((2,), interval, (x,)) == ({x: 1},)
    assert squared_distance_range(left_point, right_point, (x,)) == sp.Eq(
        sp.Symbol("t", real=True), 4
    )


def test_interval_extrema_and_level_sets_have_exact_endpoint_semantics() -> None:
    x = sp.Symbol("x", real=True)
    interval = (x >= 0) & (x <= 1)

    assert equivalent(argmax_set(x, interval, (x,)), sp.Eq(x, 1), (x,))
    assert critical_values(x**2, interval, (x,)) == (0, 1)
    assert equivalent(extrema_set(x, interval, (x,)), sp.Or(sp.Eq(x, 0), sp.Eq(x, 1)), (x,))
    assert equivalent(level_set(x**2, 1, interval), sp.Eq(x, 1), (x,))
    assert equivalent(
        sublevel_set(x, sp.Rational(1, 2), interval), (x >= 0) & (x <= sp.Rational(1, 2)), (x,)
    )
    assert equivalent(
        superlevel_set(x, sp.Rational(1, 2), interval), (x >= sp.Rational(1, 2)) & (x <= 1), (x,)
    )


def test_dimension_boundedness_density_and_path_contracts() -> None:
    x = sp.Symbol("x", real=True)
    closed = (x >= 0) & (x <= 1)
    open_interval = (x > 0) & (x < 1)
    point = sp.Eq(x, 0)

    assert is_bounded(closed, (x,)) is True
    assert is_full_dimensional(closed, (x,)) is True
    assert has_empty_interior(point, (x,)) is True
    assert is_dense_in(open_interval, closed, (x,)) is True
    assert is_path_connected(closed, (x,)) is True
    assert euler_characteristic(closed, (x,)) == 1


def test_linear_image_and_preimage_satisfy_membership_equivalence() -> None:
    x = sp.Symbol("x", real=True)
    source = (x >= 0) & (x <= 1)

    image = linear_image(source, sp.Matrix([[2]]), (x,))
    preimage = semialgebraic_preimage((2 * x,), image, (x,), target_variables=(x,))
    assert equivalent(image, (x >= 0) & (x <= 2), (x,))
    assert equivalent(preimage, source, (x,))


def test_region_regularization_and_variable_contracts() -> None:
    x = sp.Symbol("x", real=True)
    closed = SemialgebraicRegion((x >= 0) & (x <= 1), (x,))
    opened = SemialgebraicRegion((x > 0) & (x < 1), (x,))

    assert is_regular_closed_region(closed)
    assert is_regular_open_region(opened)
    assert region_variables(closed) == (x,)
    assert region_closure_interior(closed).equals_region(opened)
    assert region_interior_closure(opened).equals_region(closed)
    assert equivalent(region_regular_locus(closed), closed.formula, (x,))
    assert local_dimension(closed, {x: sp.Rational(1, 2)}) == 1
    boundary = region_boundary_result(closed)
    assert boundary.formula == sp.Or(sp.Eq(x, 0), sp.Eq(x, 1))


def test_product_measure_centroid_and_covariance_are_consistent() -> None:
    x = sp.Symbol("x", real=True)
    interval = (x >= 0) & (x <= 1)
    product = region_product(IntervalRegion(0, 1), Point((2,)))

    assert product.dimension() == 1
    assert region_measure(interval, (x,)) == 1
    assert region_centroid(interval, (x,)) == {x: sp.Rational(1, 2)}
    assert region_covariance(interval, (x,)) == sp.Matrix([[sp.Rational(1, 12)]])


def test_component_and_topology_apis_agree_on_two_intervals() -> None:
    x = sp.Symbol("x", real=True)
    region = sp.Or((x >= -2) & (x <= -1), (x >= 1) & (x <= 2))

    assert connected_component_count(region, (x,)) == 2
    samples = connected_component_samples(region, (x,))
    assert len(samples) == 2 and {sp.sign(sample[x]) for sample in samples} == {-1, 1}
    components = component_decomposition(region, (x,))
    assert components.count == 2
    summary = topology_summary(region, (x,))
    assert summary.connected_components == 2
    assert betti_number(region, 0, (x,)) == 2
    assert roadmap(region, (x,)).component_count == 2


def test_dimension_strata_triangulation_and_hardt_results_are_certified() -> None:
    x, y, a = sp.symbols("x y a", real=True)
    mixed = sp.Or((x > 0) & (x < 1) & (y > 0) & (y < 1), sp.Eq(x, 2) & sp.Eq(y, 0))
    strata = dimension_strata(mixed, (x, y))
    triangulation = triangulate_region((x >= 0) & (x <= 1), (x,))
    hardt = hardt_trivialization(sp.Eq(y**2, a), (a,), y)

    assert {item.dimension for item in strata.strata} == {0, 2}
    assert triangulation.certified and triangulation.dimension == 1
    assert hardt.certified and {item.fiber_component_count for item in hardt.strata} == {1, 2}


def test_parametric_cover_and_map_degree_certify_simple_maps() -> None:
    u = sp.Symbol("u", real=True)
    cover = bounded_parametric_cover(IntervalRegion(0, 1))
    degree = parametric_map_degree((u**2,), (u,))

    assert cover.certified_dimension() == 1
    assert degree.degree == 2 and degree.generically_finite


def test_singularity_tangent_and_thom_contracts() -> None:
    x, y = sp.symbols("x y", real=True)
    cusp = y**2 - x**3

    assert is_singular((cusp,), (0, 0), (x, y)) is True
    assert is_smooth((x**2 + y**2 - 1,), (x, y)) is True
    assert tangent_dimension((cusp,), (0, 0), (x, y)) == 2
    positive = thom_encoding(x**2 - 2, sp.sqrt(2), x)
    encodings = thom_encodings((x + 1) * (x - 2), x)
    assert positive.sign_of(x) == 1
    assert tuple(item.root for item in encodings) == (-1, 2)
