import sympy as sp

from semialg import (
    Ellipsoid,
    HalfSpace,
    HRepresentation,
    Point,
    Polytope,
    Simplex,
    affine_image,
    affine_preimage,
)


def test_affine_image_preserves_point_and_simplex_structure():
    A = ((2, 0), (0, 3))
    assert affine_image(Point((1, 2)), A, (4, 5)).coordinates == (6, 11)
    tri = Simplex(((0, 0), (1, 0), (0, 1)))
    image = affine_image(tri, A, (4, 5))
    assert isinstance(image, Simplex)
    assert set(image.vertices) == {(4, 5), (6, 5), (4, 8)}


def test_halfspace_image_and_preimage_round_trip():
    region = HalfSpace((1, 0), (2, 0))
    A = sp.Matrix(((2, 1), (0, 1)))
    image = affine_image(region, A, (3, -1))
    restored = affine_preimage(image, A, (3, -1))
    x, y = sp.symbols("x y", real=True)
    assert sp.simplify(restored.as_formula((x, y))) == sp.simplify(region.as_formula((x, y)))


def test_h_representation_image_transforms_inequalities_exactly():
    hrep = HRepresentation(((1, 0), (-1, 0), (0, 1), (0, -1)), (1, 0, 1, 0))
    image = affine_image(hrep, ((2, 0), (0, 3)), (5, 7))
    assert set(image.vertices()) == {(5, 7), (7, 7), (5, 10), (7, 10)}


def test_ellipsoid_affine_image_uses_shape_congruence():
    ellipsoid = Ellipsoid((0, 0), sp.eye(2))
    image = affine_image(ellipsoid, ((2, 0), (0, 3)))
    assert isinstance(image, Ellipsoid)
    assert image.shape_matrix == sp.ImmutableMatrix(((4, 0), (0, 9)))


def test_singular_simplex_image_degrades_to_polytope():
    tri = Simplex(((0, 0), (1, 0), (0, 1)))
    image = affine_image(tri, ((1, 0), (0, 0)))
    assert isinstance(image, Polytope)
    assert image.dimension() == 1


def test_affine_transformations_reject_nonreal_data_at_boundary() -> None:
    import pytest

    from semialg import analyze_affine_map
    from semialg.standard_regions import Point

    point = Point((sp.Integer(0),))
    message = "matrix and offset entries must be real-valued"
    with pytest.raises(ValueError, match=message):
        affine_image(point, [[sp.I]])
    with pytest.raises(ValueError, match=message):
        affine_preimage(point, [[sp.I]])
    with pytest.raises(ValueError, match=message):
        analyze_affine_map([[sp.I]])
