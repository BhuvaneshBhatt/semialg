from __future__ import annotations

import sympy as sp

from semialg import simplify_under_assumptions
from semialg.reasoning import AssumptionSimplificationResult


def test_shifted_sqrt_square_under_assumptions():
    x = sp.symbols("x", real=True)
    assert simplify_under_assumptions(sp.sqrt((x - 1) ** 2), x >= 1, [x]) == x - 1
    assert simplify_under_assumptions(sp.sqrt((x - 1) ** 2), x <= 1, [x]) == 1 - x


def test_product_sqrt_square_under_sign_assumptions():
    x, y = sp.symbols("x y", real=True)
    assert simplify_under_assumptions(sp.sqrt(x**2 * y**2), (x >= 0) & (y <= 0), [x, y]) == -x * y


def test_log_rewrites_are_domain_aware():
    x = sp.symbols("x", real=True)
    assert simplify_under_assumptions(sp.log(sp.exp(x)), True, [x]) == x
    assert simplify_under_assumptions(sp.log(x**2), x > 0, [x]) == 2 * sp.log(x)
    assert simplify_under_assumptions(sp.log(x**2), x < 0, [x]) == 2 * sp.log(-x)


def test_rational_cancellation_requires_proved_nonzero_denominator():
    x = sp.symbols("x", real=True)
    expr = (x**2 - 1) / (x - 1)
    assert simplify_under_assumptions(expr, x > 1, [x]) == x + 1
    assert simplify_under_assumptions(expr, True, [x]) == expr


def test_rational_cancellation_can_return_side_conditions():
    x = sp.symbols("x", real=True)
    expr = (x**2 - 1) / (x - 1)
    result = simplify_under_assumptions(expr, True, [x], return_conditions=True)
    assert isinstance(result, AssumptionSimplificationResult)
    assert result.expression == x + 1
    assert sp.Ne(x - 1, 0) in result.conditions
    assert "cancel_with_side_condition" in result.rewrites


def test_result_object_records_rewrites():
    x = sp.symbols("x", real=True)
    result = simplify_under_assumptions(sp.sqrt(x**2), x >= 0, [x], return_result=True)
    assert result.expression == x
    assert "sqrt_square" in result.rewrites or "abs_sign" in result.rewrites
