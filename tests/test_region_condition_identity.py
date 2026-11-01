"""Canonical region shortcuts must preserve parameter validity domains."""

import pytest
import sympy as sp

from semialg import (
    Polygon,
    Polytope,
    Simplex,
    TransformedRegion,
    affine_image,
    as_semialgebraic_region,
    is_disjoint,
    is_equal,
    is_subset,
)


def test_conditioned_simplex_relations():
    p = sp.Symbol("p", real=True)
    vertices = ((0, 0), (1, 0), (0, 1))
    left = Simplex(vertices, construction_conditions=(p > 0,))
    right = Simplex(vertices, construction_conditions=(p < 0,))
    assert left != right
    assert not is_equal(left, right)
    assert not is_subset(left, right)
    assert is_disjoint(left, right)


def test_affine_image_retains_conditions():
    p, x, y = sp.symbols("p x y", real=True)
    source = Polygon(((0, 0), (1, 0), (0, 1)), construction_conditions=(p > 0,))
    image = affine_image(source, ((2, 0), (0, 2)))
    assert isinstance(image, TransformedRegion)
    formula = as_semialgebraic_region(image, (x, y)).quantifier_free_formula()
    assert formula.subs({p: 1, x: sp.Rational(1, 2), y: sp.Rational(1, 2)}) is sp.true
    assert formula.subs({p: -1, x: sp.Rational(1, 2), y: sp.Rational(1, 2)}) is sp.false


def test_polytope_conversion_rejects_concavity():
    source = Polygon(((0, 0), (3, 0), (3, 1), (1, 1), (1, 3), (0, 3)))
    with pytest.raises(ValueError, match="concave"):
        Polytope.from_region(source)


def test_polytope_conversion_rejects_unresolved_domain():
    p = sp.Symbol("p", real=True)
    source = Simplex(((0, 0), (1, 0), (0, 1)), construction_conditions=(p > 0,))
    with pytest.raises(ValueError, match="conditions"):
        Polytope.from_region(source)


@pytest.mark.parametrize("kind", [Simplex, Polygon])
def test_copy_preserves_conditions(kind):
    p = sp.Symbol("p", real=True)
    source = kind(((0, 0), (1, 0), (0, 1)), construction_conditions=(p > 0,))
    assert kind(source) == source
    assert kind.from_region(source) == source
    copied = kind(source, construction_conditions=(p < 2,))
    assert copied.conditions == (p > 0, p < 2)
