from __future__ import annotations

import sympy as sp

from semialg import BoxRegion, SemialgebraicRegion, as_semialgebraic_region


def test_equivalent_interval_formulations_have_same_region_answers():
    x = sp.symbols("x", real=True)
    a = SemialgebraicRegion(x**2 <= 1, (x,))
    b = SemialgebraicRegion(sp.And(x >= -1, x <= 1), (x,))

    assert a.equals_region(b)
    assert a.contains((0,)) == b.contains((0,))
    assert a.contains((2,)) == b.contains((2,))
    assert a.euler_characteristic() == b.euler_characteristic() == 1
    assert sp.simplify(a.measure() - b.measure()) == 0


def test_standard_box_and_formula_region_are_semantically_identical():
    x, y = sp.symbols("x y", real=True)
    standard = as_semialgebraic_region(BoxRegion(((0, 1), (-2, 2))), (x, y))
    formula = SemialgebraicRegion(sp.And(x >= 0, x <= 1, y >= -2, y <= 2), (x, y))

    assert standard.equals_region(formula)
    assert sp.simplify(standard.measure() - formula.measure()) == 0


def test_affine_equivalent_intervals_preserve_euler_and_scaled_measure():
    x, y = sp.symbols("x y", real=True)
    source = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))
    image = source.image(2 * x + 3, variables=(y,))
    expected = SemialgebraicRegion(sp.And(y >= 3, y <= 5), (y,))

    assert image.equals_region(expected)
    assert image.euler_characteristic() == source.euler_characteristic() == 1
    assert sp.simplify(image.measure() - 2 * source.measure()) == 0
