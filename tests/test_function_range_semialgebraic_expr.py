from __future__ import annotations

import sympy as sp

from semialg import equivalent, function_range


def test_function_range_abs_expression() -> None:
    x, t = sp.symbols("x t", real=True)

    result = function_range(sp.Abs(x), True, [x], value_symbol=t, return_result=True)

    assert equivalent(result.range_condition, t >= 0, [t])
    assert result.lower_bound == 0
    assert result.upper_bound == sp.oo
    assert result.lower_bound_attained is True
    assert result.is_interval is True


def test_function_range_sqrt_expression() -> None:
    x, t = sp.symbols("x t", real=True)

    result = function_range(sp.sqrt(1 - x**2), True, [x], value_symbol=t, return_result=True)

    assert equivalent(result.range_condition, sp.And(t >= 0, t <= 1), [t])
    assert result.lower_bound == 0
    assert result.upper_bound == 1
    assert result.lower_bound_attained is True
    assert result.upper_bound_attained is True


def test_function_range_max_and_min_expressions() -> None:
    x, t = sp.symbols("x t", real=True)

    max_result = function_range(
        sp.Max(x, 0), sp.And(x >= -1, x <= 2), [x], value_symbol=t, return_result=True
    )
    min_result = function_range(
        sp.Min(x, 1), sp.And(x >= 0, x <= 3), [x], value_symbol=t, return_result=True
    )

    assert equivalent(max_result.range_condition, sp.And(t >= 0, t <= 2), [t])
    assert equivalent(min_result.range_condition, sp.And(t >= 0, t <= 1), [t])


def test_function_range_piecewise_expression() -> None:
    x, t = sp.symbols("x t", real=True)
    expr = sp.Piecewise((x**2, x < 0), (x, True))

    result = function_range(expr, sp.And(x >= -2, x <= 3), [x], value_symbol=t, return_result=True)

    assert equivalent(result.range_condition, sp.And(t >= 0, t <= 4), [t])
    assert result.lower_bound == 0
    assert result.upper_bound == 4


def test_function_range_semialgebraic_subexpression_inside_arithmetic() -> None:
    x, t = sp.symbols("x t", real=True)

    result = function_range(sp.Abs(x) + 1, True, [x], value_symbol=t, return_result=True)

    assert equivalent(result.range_condition, t >= 1, [t])
    assert result.lower_bound == 1
    assert result.upper_bound == sp.oo
