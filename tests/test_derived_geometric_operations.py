import sympy as sp

from semialg import (
    affine_transform,
    argmax_set,
    argmin_set,
    centroid,
    closest_points,
    connected_components,
    contains_point,
    coordinate_range,
    covariance_matrix,
    diameter,
    extrema_set,
    has_empty_interior,
    inertia_tensor,
    intersects,
    is_bounded,
    is_closed,
    is_compact,
    is_connected,
    is_disjoint,
    is_empty,
    is_equal,
    is_full_dimensional,
    is_open,
    is_subset,
    level_set,
    minkowski_sum,
    moment_matrix,
    nearest_point,
    scale,
    sublevel_set,
    superlevel_set,
    support_function,
    translate,
    width,
)


def test_global_extrema_sets_retain_positive_dimensional_optimizers():
    x, y = sp.symbols("x y", real=True)
    box = sp.And(x >= -1, x <= 1, y >= -2, y <= 2)
    minimum = argmin_set(x**2, box, [x, y])
    maximum = argmax_set(x**2, box, [x, y])
    both = extrema_set(x**2, box, [x, y])
    assert contains_point(minimum, {x: 0, y: 1}, [x, y])
    assert not contains_point(minimum, {x: 1, y: 1}, [x, y])
    assert contains_point(maximum, {x: -1, y: 0}, [x, y])
    assert contains_point(maximum, {x: 1, y: 0}, [x, y])
    assert contains_point(both, {x: 0, y: 0}, [x, y])
    assert contains_point(both, {x: 1, y: 0}, [x, y])


def test_level_set_constructors_preserve_region_constraint():
    x = sp.symbols("x", real=True)
    region = x >= 0
    assert contains_point(level_set(x**2, 1, region), {x: 1}, [x])
    assert not contains_point(level_set(x**2, 1, region), {x: -1}, [x])
    assert contains_point(sublevel_set(x**2, 1, region), {x: sp.Rational(1, 2)}, [x])
    assert not contains_point(sublevel_set(x**2, 1, region), {x: 2}, [x])
    assert contains_point(superlevel_set(x**2, 1, region), {x: 2}, [x])


def test_geometric_property_and_set_relation_api():
    x = sp.symbols("x", real=True)
    closed = sp.And(x >= 0, x <= 1)
    opened = sp.And(x > 0, x < 1)
    assert is_bounded(closed, [x])
    assert is_compact(closed, [x])
    assert is_closed(closed, [x])
    assert is_open(opened, [x])
    assert not is_empty(closed, [x])
    assert is_subset(opened, closed, [x])
    assert not is_equal(opened, closed, [x])
    assert is_disjoint(x < 0, x > 0, [x])
    assert intersects(closed, x > sp.Rational(1, 2), [x])
    assert is_connected(closed, [x])
    assert is_full_dimensional(closed, [x])
    assert has_empty_interior(sp.Eq(x, 0), [x])


def test_coordinate_range_diameter_nearest_and_support_operations():
    x = sp.symbols("x", real=True)
    interval = sp.And(x >= 0, x <= 2)
    t = sp.Symbol("t", real=True)
    rng = coordinate_range(interval, x, [x])
    assert sp.simplify(rng.subs(t, 1)) is sp.true
    assert sp.simplify(rng.subs(t, 3)) is sp.false
    assert diameter(interval, [x]) == 2
    assert support_function(interval, (1,), [x]) == 2
    assert width(interval, (1,), [x]) == 2
    nearest = nearest_point((3,), interval, [x])
    assert nearest and nearest[0][x] == 2
    pairs = closest_points(x <= 0, x >= 2, [x])
    assert pairs and pairs[0][0][x] == 0 and pairs[0][1][x] == 2


def test_moment_centroid_covariance_and_inertia_conveniences():
    x, y = sp.symbols("x y", real=True)
    square = sp.And(x >= -1, x <= 1, y >= -1, y <= 1)
    assert centroid(square, [x, y]) == {x: 0, y: 0}
    assert moment_matrix(square, [x, y]) == sp.diag(sp.Rational(1, 3), sp.Rational(1, 3))
    assert covariance_matrix(square, [x, y]) == sp.diag(sp.Rational(1, 3), sp.Rational(1, 3))
    assert inertia_tensor(square, [x, y]) == sp.diag(sp.Rational(4, 3), sp.Rational(4, 3))


def test_affine_transform_translate_scale_and_minkowski_sum():
    x = sp.symbols("x", real=True)
    interval = sp.And(x >= 0, x <= 1)
    assert is_equal(translate(interval, (2,), [x]), sp.And(x >= 2, x <= 3), [x])
    assert is_equal(scale(interval, 2, [x]), sp.And(x >= 0, x <= 2), [x])
    affine = affine_transform(interval, [[2]], [1], [x])
    assert is_equal(affine, sp.And(x >= 1, x <= 3), [x])
    summed = minkowski_sum(interval, sp.And(x >= 2, x <= 4), [x])
    assert is_equal(summed, sp.And(x >= 2, x <= 5), [x])


def test_dense_linear_image_and_distance_set_operations():
    x = sp.symbols("x", real=True)
    from semialg import distance_set, is_dense_in, linear_image, squared_distance_range

    assert is_dense_in(sp.Ne(x, 0), sp.true, [x])
    image = linear_image(sp.And(x >= 0, x <= 1), [[3]], [x])
    assert is_equal(image, sp.And(x >= 0, x <= 3), [x])
    squared = squared_distance_range(sp.And(x >= 0, x <= 2), variables=[x])
    t = next(sym for sym in squared.free_symbols if sym.name == "t")
    assert sp.simplify(squared.subs(t, 4)) is sp.true
    distances = distance_set(sp.And(x >= 0, x <= 2), variables=[x])
    d = next(sym for sym in distances.free_symbols if sym.name == "d")
    assert sp.simplify(distances.subs(d, 2)) is sp.true
    assert sp.simplify(distances.subs(d, -1)) is sp.false


def test_connected_components_public_api_returns_exact_components():
    x = sp.Symbol("x", real=True)
    components = connected_components((x <= -1) | (x >= 1), [x])

    assert len(components) == 2
    assert any(sp.simplify(component.subs(x, -2)) is sp.true for component in components)
    assert any(sp.simplify(component.subs(x, 2)) is sp.true for component in components)


def test_connected_components_separates_strict_intervals_at_an_excluded_point():
    x = sp.Symbol("x", real=True)

    components = connected_components((x < 0) | (x > 0), [x])

    assert components == (x < 0, x > 0)
