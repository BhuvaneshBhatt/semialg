import sympy as sp

from semialg import (
    Ball,
    Ellipsoid,
    Polytope,
    SemialgebraicRegion,
    Simplex,
    TransformedRegion,
    region_image,
    region_preimage,
)


def test_region_image_detects_affine_map_and_preserves_simplex():
    x, y = sp.symbols("x y", real=True)
    triangle = Simplex(((0, 0), (1, 0), (0, 1)))
    image = region_image(triangle, (2 * x + y + 3, y - 1), (x, y))
    assert isinstance(image, Simplex)
    assert set(image.vertices) == {(3, -1), (5, -1), (4, 0)}


def test_region_image_keeps_nonlinear_map_structural_until_lowering():
    x, y = sp.symbols("x y", real=True)
    disk = Ball((0, 0), 1)
    image = region_image(disk, (x, y**2), (x, y))
    assert isinstance(image, TransformedRegion)
    assert image.ambient_dimension() == 2


def test_region_preimage_detects_invertible_affine_ellipsoid_map():
    x, y = sp.symbols("x y", real=True)
    target = Ellipsoid((0, 0), sp.eye(2))
    preimage = region_preimage(target, (2 * x, 3 * y), (x, y))
    assert isinstance(preimage, Ellipsoid)
    assert preimage.shape_matrix == sp.diag(sp.Rational(1, 4), sp.Rational(1, 9))


def test_region_preimage_of_nonlinear_map_returns_semialgebraic_region():
    x, y = sp.symbols("x y", real=True)
    target = Ball((0, 0), 1)
    preimage = region_preimage(target, (x, y**2), (x, y))
    assert isinstance(preimage, SemialgebraicRegion)
    assert sp.simplify(preimage.formula) == (x**2 + y**4 <= 1)


def test_region_image_affine_polytope_preserves_combinatorics():
    x, y = sp.symbols("x y", real=True)
    square = Polytope(((0, 0), (1, 0), (1, 1), (0, 1)))
    image = region_image(square, (x + y, 2 * y), (x, y))
    assert isinstance(image, Polytope)
    assert square.is_combinatorially_equivalent(image)
