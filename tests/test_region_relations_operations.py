import sympy as sp

from semialg import (
    Ball,
    BoxRegion,
    IntervalRegion,
    Point,
    SemialgebraicRegion,
    is_disjoint,
    is_equal,
    is_subset,
    minkowski_sum,
    region_image,
    region_intersection,
    region_preimage,
    region_product,
)


def test_interval_relations_are_exact():
    inner = IntervalRegion(0, 1)
    outer = IntervalRegion(-1, 2)
    touching = IntervalRegion(1, 3)
    separate = IntervalRegion(2, 3)

    assert is_subset(inner, outer)
    assert not is_disjoint(inner, touching)
    assert is_disjoint(inner, separate)
    assert is_equal(inner, IntervalRegion(0, 1))


def test_box_relations_are_coordinatewise():
    inner = BoxRegion(((0, 1), (0, 2)))
    outer = BoxRegion(((-1, 3), (-2, 4)))
    separate = BoxRegion(((2, 3), (0, 2)))

    assert inner.subset_of(outer)
    assert inner.disjoint_from(separate)
    assert not inner.equals_region(outer)


def test_ball_relations_use_euclidean_geometry():
    inner = Ball((1, 0), 1)
    outer = Ball((0, 0), 3)
    touching = Ball((4, 0), 1)
    separate = Ball((5, 0), 1)

    assert inner.subset_of(outer)
    assert not outer.disjoint_from(touching)
    assert outer.disjoint_from(separate)


def test_structural_intersection_preserves_region_type():
    result = region_intersection(IntervalRegion(0, 3), IntervalRegion(1, 2))
    assert result == IntervalRegion(1, 2)

    box = region_intersection(
        BoxRegion(((0, 3), (-2, 2))),
        BoxRegion(((1, 4), (-1, 1))),
    )
    assert box == BoxRegion(((1, 3), (-1, 1)))


def test_cartesian_product_preserves_boxes_and_points():
    box = region_product(IntervalRegion(0, 1), IntervalRegion(2, 4))
    assert box == BoxRegion(((0, 1), (2, 4)))

    point = region_product(Point((1,)), Point((2, 3)))
    assert point == Point((1, 2, 3))


def test_minkowski_sum_preserves_common_canonical_regions():
    assert minkowski_sum(IntervalRegion(0, 1), IntervalRegion(2, 4)) == IntervalRegion(2, 5)
    assert minkowski_sum(BoxRegion(((0, 1), (1, 2))), BoxRegion(((-1, 2), (3, 5)))) == BoxRegion(
        ((-1, 3), (4, 7))
    )
    assert minkowski_sum(Ball((0, 1), 2), Ball((3, -1), 4)) == Ball((3, 0), 6)


def test_point_minkowski_sum_is_translation():
    result = Point((2, -1)).minkowski_sum(Ball((0, 0), 3))
    assert result == Ball((2, -1), 3)


def test_general_intersection_falls_back_to_symbolic_region():
    x, y = sp.symbols("x y", real=True)
    disk = Ball((0, 0), 2)
    half_plane = SemialgebraicRegion(x >= 0, (x, y))

    result = region_intersection(disk, half_plane)

    assert isinstance(result, SemialgebraicRegion)
    assert result.contains((1, 0))
    assert not result.contains((-1, 0))


def test_polynomial_preimage_is_exact():
    x = sp.symbols("x", real=True)
    target = IntervalRegion(0, 1)

    preimage = region_preimage(target, (x**2,), (x,))

    assert isinstance(preimage, SemialgebraicRegion)
    assert preimage.contains((-1,))
    assert preimage.contains((sp.Rational(1, 2),))
    assert not preimage.contains((2,))


def test_polynomial_image_lowers_to_semialgebraic_formula():
    x, y = sp.symbols("x y", real=True)
    source = IntervalRegion(-1, 1)

    image = region_image(source, (x**2,), (x,))
    symbolic = image.as_semialgebraic_region((y,))

    assert symbolic.contains((0,))
    assert symbolic.contains((1,))
    assert not symbolic.contains((2,))


def test_empty_open_singleton_relations():
    empty = IntervalRegion(0, 0, lower_closed=False, upper_closed=False)
    nonempty = IntervalRegion(1, 2)

    assert empty.subset_of(nonempty)
    assert empty.disjoint_from(nonempty)
    assert region_intersection(empty, nonempty) == empty


def test_open_interval_product_preserves_endpoint_semantics():
    x, y = sp.symbols("x y", real=True)
    result = region_product(
        IntervalRegion(0, 1, lower_closed=False),
        IntervalRegion(2, 3, upper_closed=False),
        variables=(x, y),
    )

    assert isinstance(result, SemialgebraicRegion)
    assert result.contains((sp.Rational(1, 2), sp.Rational(5, 2)))
    assert not result.contains((0, sp.Rational(5, 2)))
    assert not result.contains((sp.Rational(1, 2), 3))
