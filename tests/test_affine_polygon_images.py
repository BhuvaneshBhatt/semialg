"""Affine embeddings retain concavities instead of taking a convex hull."""

import sympy as sp

from semialg import Polygon, Polytope, TransformedRegion, affine_image, as_semialgebraic_region


def test_concave_embedding():
    polygon = Polygon(((0, 0), (3, 0), (3, 1), (1, 1), (1, 3), (0, 3)))
    image = affine_image(polygon, ((1, 0), (0, 1), (0, 0)))
    assert isinstance(image, TransformedRegion)
    x, y, z = sp.symbols("x y z", real=True)
    region = as_semialgebraic_region(image, (x, y, z))
    assert region.contains((sp.Rational(1, 2), 2, 0))
    assert not region.contains((sp.Rational(3, 2), sp.Rational(3, 2), 0))
    assert not region.contains((sp.Rational(1, 2), 2, 1))


def test_convex_embedding():
    polygon = Polygon(((0, 0), (1, 0), (1, 1), (0, 1)))
    image = affine_image(polygon, ((1, 0), (0, 1), (0, 0)))
    assert isinstance(image, Polytope)
    assert image.dimension() == 2
    assert set(image.vertices) == {(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)}
