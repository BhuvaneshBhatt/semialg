import pytest
import sympy as sp

from semialg import (
    canonicalize_one_dimensional_formula,
    function_domain,
    is_real_valued,
    normalize_domain_sensitive_constraints,
    solve_semialgebraic,
)

pytestmark = pytest.mark.slow


def test_function_domain_and_real_valuedness_helpers():
    x = sp.symbols("x", real=True)
    assert function_domain(sp.sqrt(x - 1), [x]) == (x - 1 >= 0)
    assert is_real_valued(sp.sqrt(x - 1), [x], assumptions=x >= 1) is True
    assert is_real_valued(sp.sqrt(x - 1), [x], assumptions=x < 1) is False


def test_domain_sensitive_rational_solving_rewrites_denominator_safely():
    x = sp.symbols("x", real=True)
    normalized = normalize_domain_sensitive_constraints(1 / (x - 1) > 0, [x]).formula
    assert normalized == (x > 1)
    result = solve_semialgebraic([1 / (x - 1) > 0], [x], count=0, output="reduced_formula")
    assert result == (x > 1)


def test_domain_sensitive_sqrt_solving():
    x = sp.symbols("x", real=True)
    result = solve_semialgebraic([sp.sqrt(x - 1) <= 2], [x], count=0, output="reduced_formula")
    assert result == sp.And(x >= 1, x <= 5)


def test_canonical_one_dimensional_formula_output():
    x = sp.symbols("x", real=True)
    result = canonicalize_one_dimensional_formula(sp.And(x**2 <= 1, sp.Ne(x, 0)), x)
    assert result == sp.Or(sp.And(x >= -1, x < 0), sp.And(x > 0, x <= 1))


def test_method_and_variable_order_controls_are_recorded():
    x, y = sp.symbols("x y", real=True)
    result = solve_semialgebraic(
        [x >= 0, y >= 0], [x, y], variable_order=[y, x], method="cad", count=0
    )
    assert tuple(result.variables) == (y, x)
    assert result.diagnostics["requested_method"] == "cad"
    assert result.diagnostics["variable_order"] == ("y", "x")
