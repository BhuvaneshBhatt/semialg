from __future__ import annotations

import sympy as sp

from semialg import classify_real_roots, function_range


def test_root_classification_result_has_stable_value_equality_and_repr() -> None:
    x, a = sp.symbols("x a", real=True)
    first = classify_real_roots(x**2 + a, x, parameters=[a])
    second = classify_real_roots(x**2 + a, x, parameters=[a])

    assert first == second
    assert repr(first) == repr(second)


def test_function_range_result_repr_contains_exact_formula() -> None:
    x, t = sp.symbols("x t", real=True)
    result = function_range(x**2, x >= 0, [x], value_symbol=t, return_result=True)

    rendered = repr(result)
    assert "FunctionRangeResult" in rendered
    assert "formula=" in rendered
