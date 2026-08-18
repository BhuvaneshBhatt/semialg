from __future__ import annotations

import sympy as sp

from semialg import MeasureResult, semialgebraic_measure


def test_measure_1d_interval_from_polynomial_inequality() -> None:
    x = sp.symbols("x", real=True)
    assert semialgebraic_measure(x**2 <= 1, [x]) == 2


def test_measure_1d_boolean_combination_with_bounds() -> None:
    x = sp.symbols("x", real=True)
    value = semialgebraic_measure((x < -1) | (x > 1), [x], bounds=[(x, -3, 3)])
    assert value == 4


def test_measure_1d_zero_measure_equality() -> None:
    x = sp.symbols("x", real=True)
    assert semialgebraic_measure(sp.Eq(x, 0), [x]) == 0


def test_measure_disk_and_annulus() -> None:
    x, y = sp.symbols("x y", real=True)
    assert semialgebraic_measure(x**2 + y**2 <= 1, [x, y]) == sp.pi
    assert semialgebraic_measure(sp.And(x**2 + y**2 <= 4, x**2 + y**2 >= 1), [x, y]) == 3 * sp.pi


def test_measure_triangle_by_vertical_slices() -> None:
    x, y = sp.symbols("x y", real=True)
    region = sp.And(x >= 0, y >= 0, x + y <= 1)
    assert semialgebraic_measure(region, [x, y]) == sp.Rational(1, 2)


def test_measure_parabola_cap_by_vertical_slices() -> None:
    x, y = sp.symbols("x y", real=True)
    region = sp.And(y >= x**2, y <= 1)
    assert semialgebraic_measure(region, [x, y]) == sp.Rational(4, 3)


def test_measure_result_object_records_method() -> None:
    x = sp.symbols("x", real=True)
    result = semialgebraic_measure(x**2 <= 1, [x], return_result=True)
    assert isinstance(result, MeasureResult)
    assert result.value == 2
    assert result.method == "one_dimensional_cell_sampling"
