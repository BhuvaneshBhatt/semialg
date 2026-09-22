import sympy as sp

from semialg import Ball, Box, Point, Simplex, Sphere


def _mapped_point(matrix, offset, point):
    value = sp.Matrix(matrix) * sp.Matrix(point) + sp.Matrix(offset)
    return tuple(sp.simplify(item) for item in value)


def test_full_dimensional_measure_scales_by_determinant():
    region = Simplex(((0, 0), (2, 0), (0, 3)))
    matrix = sp.Matrix(((2, 1), (0, 3)))
    offset = (5, -2)
    image = region.transform(matrix, offset)
    expected = sp.Abs(matrix.det()) * region.measure()
    assert sp.simplify(image.measure() - expected) == 0


def test_centroid_is_affine_equivariant():
    region = Box(((-1, 3), (2, 6)))
    matrix = sp.Matrix(((2, -1), (1, 1)))
    offset = (4, -3)
    image = region.transform(matrix, offset)
    expected = _mapped_point(matrix, offset, region.centroid())
    assert image.centroid() == expected


def test_segment_measure_uses_intrinsic_scale():
    region = Simplex(((0, 0), (1, 0)))
    matrix = sp.Matrix(((3, 0), (4, 2)))
    image = region.transform(matrix, (0, 0))
    assert sp.simplify(image.measure() - 5) == 0
    assert image.measure(measure_dimension="ambient") == 0


def test_similarity_preserves_ball_family_and_measure_law():
    region = Ball((1, -2), 3)
    matrix = sp.Matrix(((0, -2), (2, 0)))
    image = region.transform(matrix, (5, 7))
    assert isinstance(image, Ball)
    assert image.center == (9, 9)
    assert sp.simplify(image.measure() - 4 * region.measure()) == 0


def test_anisotropic_ball_image_is_ellipsoid_with_correct_volume():
    region = Ball((0, 0), 1)
    matrix = sp.diag(2, 3)
    image = region.transform(matrix, (1, -1))
    assert image.center == (1, -1)
    assert image.shape_matrix == sp.diag(4, 9)
    assert sp.simplify(image.measure() - 6 * sp.pi) == 0


def test_sphere_similarity_scales_intrinsic_measure():
    region = Sphere((0, 0, 0), 2)
    matrix = 3 * sp.eye(3)
    image = region.transform(matrix, (1, 2, 3))
    assert sp.simplify(image.measure() - 9 * region.measure()) == 0


def test_point_transform_preserves_counting_measure():
    region = Point((2, -1))
    image = region.transform(((2, 0), (0, 3)), (4, 5))
    assert image.coordinates == (8, 2)
    assert image.measure() == region.measure() == 1


def test_symbolic_rank_is_not_assumed_invertible():
    a = sp.Symbol("a", real=True)
    region = Ball((0, 0), 1)
    image = region.transform(sp.diag(a, 1), (0, 0))
    from semialg import TransformedRegion

    assert isinstance(image, TransformedRegion)


def test_nonzero_symbolic_scale_preserves_canonical_image():
    a = sp.Symbol("a", real=True, nonzero=True)
    region = Ball((0, 0), 1)
    image = region.transform(sp.diag(a, 1), (0, 0))
    assert image.shape_matrix == sp.diag(a**2, 1)
