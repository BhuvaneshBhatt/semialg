from __future__ import annotations

import sympy as sp

from semialg import FunctionRangeResult, function_range


def test_function_range_identity_over_disconnected_domain_is_disconnected() -> None:
    x, t = sp.symbols("x t", real=True)
    result = function_range(x, sp.Or(x <= -1, x >= 1), [x], value_symbol=t, return_result=True)
    assert isinstance(result, FunctionRangeResult)
    assert result.method == "qe_image_solved_graph"
    assert sp.simplify(result.formula) == (sp.Abs(t) >= 1)


def test_function_range_affine_open_interval_preserves_strictness() -> None:
    x, t = sp.symbols("x t", real=True)
    formula = function_range(2 * x + 1, sp.And(x > 0, x < 1), [x], value_symbol=t)
    assert sp.simplify(formula) == sp.And(t > 1, t < 3)


def test_function_range_univariate_rational_positive_ray() -> None:
    x, t = sp.symbols("x t", real=True)
    formula = function_range(1 / x, x > 0, [x], value_symbol=t)
    assert sp.simplify(formula) == (t > 0)


def test_function_range_falls_back_to_bounds_for_quadratic_interval() -> None:
    x, t = sp.symbols("x t", real=True)
    result = function_range(x**2, [x >= -2, x <= 3], [x], value_symbol=t, return_result=True)
    assert isinstance(result, FunctionRangeResult)
    assert result.formula == sp.And(t >= 0, t <= 9)
    assert "optimization_bounds" in result.method
