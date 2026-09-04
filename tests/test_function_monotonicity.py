import sympy as sp

from semialg import function_convexity, function_monotonicity
from semialg.conditional import ParameterStratifiedResult
from semialg.function_analysis import FunctionMonotonicityResult


def test_polynomial_monotonicity_uses_derivative_sign_and_strict_zero_set():
    x = sp.symbols("x", real=True)
    assert function_monotonicity(x, x) == "strictly_increasing"
    assert function_monotonicity(-x, x) == "strictly_decreasing"
    assert function_monotonicity(x**3, x) == "strictly_increasing"
    assert function_monotonicity(x**2, x) == "nonmonotonic"
    assert function_monotonicity(5, x) == "constant"


def test_nonsmooth_abs_is_nonmonotonic_globally_and_monotone_on_halfline():
    x = sp.symbols("x", real=True)
    assert function_monotonicity(sp.Abs(x), x) == "nonmonotonic"
    assert function_monotonicity(sp.Abs(x), x, domain=x >= 0) == "strictly_increasing"


def test_structured_monotonicity_result_contains_certificates():
    x = sp.symbols("x", real=True)
    result = function_monotonicity(x**3, x, return_result=True)
    assert isinstance(result, FunctionMonotonicityResult)
    assert result.certified
    assert result.classification == "strictly_increasing"
    assert result.increasing is True
    assert result.strictly_increasing is True


def test_parameterized_linear_monotonicity_is_conditional_automatically():
    x, a = sp.symbols("x a", real=True)
    result = function_monotonicity(a * x, x)
    assert isinstance(result, ParameterStratifiedResult)
    assert result.select({a: 2}) == "strictly_increasing"
    assert result.select({a: 0}) == "constant"
    assert result.select({a: -2}) == "strictly_decreasing"


def test_explicit_parameter_argument_matches_automatic_parameter_detection():
    x, a = sp.symbols("x a", real=True)
    result = function_monotonicity(a * x, x, parameters=[a])
    assert result.select({a: 1}) == "strictly_increasing"
    assert result.select({a: -1}) == "strictly_decreasing"


def test_function_convexity_automatically_returns_conditions_for_free_parameters():
    x, a = sp.symbols("x a", real=True)
    result = function_convexity(a * x**2, [x])
    assert isinstance(result, ParameterStratifiedResult)
    assert result.select({a: 1}) == "strongly_convex"
    assert result.select({a: 0}) == "affine"
    assert result.select({a: -1}) == "strongly_concave"
