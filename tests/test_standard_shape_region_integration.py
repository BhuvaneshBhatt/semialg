from __future__ import annotations

import pytest
import sympy as sp

from semialg import integrate_over_region, reduce_region_integral, semialgebraic_measure

pytestmark = pytest.mark.slow


def test_axis_aligned_box_reduction_and_integration() -> None:
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x >= 0, x <= 2, y >= 1, y <= 3)

    reduced = reduce_region_integral(x + y, condition, [x, y])

    assert reduced.method == "axis_aligned_box_iterated_integral"
    assert len(reduced.pieces) == 1
    assert integrate_over_region(x + y, condition, [x, y]) == 12
    assert semialgebraic_measure(condition, [x, y]) == 4


def test_unit_simplex_standard_shape() -> None:
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x >= 0, y >= 0, x + y <= 1)

    reduced = reduce_region_integral(x * y, condition, [x, y])

    assert reduced.method == "unit_simplex_iterated_integral"
    assert integrate_over_region(x * y, condition, [x, y]) == sp.Rational(1, 24)
    assert semialgebraic_measure(condition, [x, y]) == sp.Rational(1, 2)


def test_axis_aligned_ellipse_standard_shape() -> None:
    x, y = sp.symbols("x y", real=True)
    condition = (x - 1) ** 2 / 4 + (y + 2) ** 2 / 9 <= 1

    reduced = reduce_region_integral(1, condition, [x, y])

    assert reduced.method == "axis_aligned_ellipse_affine_unit_disk"
    assert integrate_over_region(1, condition, [x, y]) == 6 * sp.pi
    assert integrate_over_region(x, condition, [x, y]) == 6 * sp.pi
    assert semialgebraic_measure(condition, [x, y]) == 6 * sp.pi


def test_higher_dimensional_box_standard_shape() -> None:
    x, y, z = sp.symbols("x y z", real=True)
    condition = sp.And(x >= 0, x <= 1, y >= 0, y <= 1, z >= 0, z <= 1)

    reduced = reduce_region_integral(x + y + z, condition, [x, y, z])

    assert reduced.method == "axis_aligned_box_iterated_integral"
    assert integrate_over_region(x + y + z, condition, [x, y, z]) == sp.Rational(3, 2)
    assert semialgebraic_measure(condition, [x, y, z]) == 1
