import sympy as sp

from semialg import (
    AffineHalfSpace,
    Ball,
    Ellipsoid,
    EllipsoidBoundary,
    HRepresentation,
    Polygon,
    SemialgebraicRegion,
    Sphere,
    TransformedRegion,
    affine_image,
    affine_preimage,
)


def test_geometry_transform_method_routes_to_affine_image():
    image = Ball((0, 0), 2).transform(((0, -3), (3, 0)), (1, 2))
    assert isinstance(image, Ball)
    assert image.center == (1, 2)
    assert image.radius == 6


def test_anisotropic_ball_image_is_ellipsoid():
    image = affine_image(Ball((1, 2), 3), ((2, 0), (0, 4)), (5, -1))
    assert isinstance(image, Ellipsoid)
    assert image.center == (7, 7)
    assert image.shape_matrix == sp.diag(36, 144)


def test_anisotropic_sphere_image_is_ellipsoid_boundary():
    image = affine_image(Sphere((0, 0), 1), ((2, 0), (0, 3)))
    assert isinstance(image, EllipsoidBoundary)
    x, y = sp.symbols("x y", real=True)
    assert sp.simplify(image.as_formula((x, y))) == sp.Eq(x**2 / 4 + y**2 / 9, 1)


def test_rank_changing_ball_image_remains_exact_transformed_region():
    image = affine_image(Ball((0, 0), 1), ((1, 0),))
    assert isinstance(image, TransformedRegion)
    assert image.ambient_dimension() == 1


def test_affine_halfspace_preserves_intrinsic_structure_under_embedding():
    region = AffineHalfSpace((0, 0), ((1, 0),), (0, 1))
    image = affine_image(region, ((1, 0), (0, 1), (1, 1)), (2, 3, 4))
    assert isinstance(image, AffineHalfSpace)
    assert image.dimension() == 2
    assert image.ambient_dimension() == 3


def test_polygon_invertible_affine_image_stays_polygon():
    polygon = Polygon(((0, 0), (2, 0), (2, 1), (0, 1)))
    image = affine_image(polygon, ((1, 1), (0, 2)), (3, 4))
    assert isinstance(image, Polygon)
    assert image.vertices == ((3, 4), (5, 4), (6, 6), (4, 6))


def test_rectangular_affine_preimage_is_exact_semialgebraic_substitution():
    target = Ball((0, 0), 1)
    preimage = affine_preimage(target, ((1, 0, 0), (0, 1, 0)))
    assert isinstance(preimage, SemialgebraicRegion)
    u0, u1, u2 = preimage.variables
    assert sp.simplify(preimage.formula) == (u0**2 + u1**2 <= 1)
    assert u2 not in preimage.formula.free_symbols


def test_singular_affine_preimage_is_exact_semialgebraic_substitution():
    target = Ball((0, 0), 1)
    preimage = affine_preimage(target, ((1, 0), (0, 0)))
    assert isinstance(preimage, SemialgebraicRegion)
    u0, _ = preimage.variables
    assert sp.simplify(preimage.formula) == (u0**2 <= 1)


def test_hrepresentation_rectangular_preimage_is_supported():
    hrep = HRepresentation(((1,), (-1,)), (1, 0))
    preimage = affine_preimage(hrep, ((1, 0),))
    assert isinstance(preimage, SemialgebraicRegion)
    u0, _ = preimage.variables
    assert sp.simplify(preimage.formula) == sp.And(u0 >= 0, u0 <= 1)


def test_zero_radius_ball_and_sphere_images_are_points():
    from semialg import Point

    matrix = ((2, 1), (0, 3))
    assert isinstance(affine_image(Ball((1, 2), 0), matrix), Point)
    assert isinstance(affine_image(Sphere((1, 2), 0), matrix), Point)
