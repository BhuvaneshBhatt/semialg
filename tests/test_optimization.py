from __future__ import annotations

import pytest
import sympy as sp

from semialg import (
    FunctionRangeResult,
    OptimizationResult,
    function_range,
    semialgebraic_maximize,
    semialgebraic_minimize,
)

pytestmark = pytest.mark.slow


def test_minimize_univariate_closed_endpoint() -> None:
    x = sp.symbols("x", real=True)
    result = semialgebraic_minimize(x**2, [x >= 2], [x])
    assert isinstance(result, OptimizationResult)
    assert result.value == 4
    assert result.attained is True
    assert result.point == {x: 2}


def test_minimize_univariate_open_endpoint_reports_infimum() -> None:
    x = sp.symbols("x", real=True)
    result = semialgebraic_minimize(x, [x > 0], [x])
    assert isinstance(result, OptimizationResult)
    assert result.value == 0
    assert result.attained is False
    assert result.points == ()


def test_maximize_triangle_product() -> None:
    x, y = sp.symbols("x y", real=True)
    region = [x >= 0, y >= 0, x + y <= 1]
    result = semialgebraic_maximize(x * y, region, [x, y])
    assert isinstance(result, OptimizationResult)
    assert result.value == sp.Rational(1, 4)
    assert result.attained is True
    assert result.point == {x: sp.Rational(1, 2), y: sp.Rational(1, 2)}


def test_minimize_distance_to_origin_over_halfplane() -> None:
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_minimize(x**2 + y**2, [x + y >= 1], [x, y])
    assert result.value == sp.Rational(1, 2)
    assert result.attained is True
    assert result.point == {x: sp.Rational(1, 2), y: sp.Rational(1, 2)}


def test_function_range_univariate_polynomial_on_interval() -> None:
    x, t = sp.symbols("x t", real=True)
    formula = function_range(x**2, [x >= -2, x <= 3], [x], value_symbol=t)
    assert formula == sp.And(t >= 0, t <= 9)


def test_function_range_open_interval_uses_strict_endpoint() -> None:
    x, t = sp.symbols("x t", real=True)
    formula = function_range(x, [x > 0, x < 1], [x], value_symbol=t)
    assert formula == sp.And(t > 0, t < 1)


def test_function_range_result_object() -> None:
    x = sp.symbols("x", real=True)
    result = function_range(x**2, [x >= -1, x <= 1], [x], return_result=True)
    assert isinstance(result, FunctionRangeResult)
    assert result.infimum == 0
    assert result.supremum == 1
    assert result.minimum_attained is True
    assert result.maximum_attained is True
