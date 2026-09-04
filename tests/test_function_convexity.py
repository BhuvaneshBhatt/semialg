import sympy as sp

from semialg import function_convexity
from semialg.conditional import ParameterStratifiedResult
from semialg.function_analysis import FunctionConvexityResult


def test_polynomial_convex_concave_affine_and_neither():
    x = sp.symbols("x", real=True)
    assert function_convexity(x**2, (x,)) == "strongly_convex"
    assert function_convexity(-(x**2), (x,)) == "strongly_concave"
    assert function_convexity(3 * x + 2, (x,)) == "affine"
    assert function_convexity(x**3, (x,)) == "neither"


def test_certificate_uses_matrix_backend_for_polynomial():
    x, y = sp.symbols("x y", real=True)
    result = function_convexity(x**2 + y**2, (x, y), return_result=True)
    assert isinstance(result, FunctionConvexityResult)
    assert result.classification == "strongly_convex"
    assert result.convex is True
    assert result.concave is False
    assert result.convexity_certificate is not None
    assert "hessian_matrix_definiteness" in result.method


def test_lower_dimensional_affine_reduction_prevents_false_hessian_rejection():
    x, y = sp.symbols("x y", real=True)
    domain = sp.And(sp.Eq(y, 0), -1 <= x, x <= 1)
    result = function_convexity(-(y**2), (x, y), domain=domain, return_result=True)
    assert result.classification == "affine"
    assert result.details["relative_strict_feasibility"].feasible is True
    assert result.details["presolve"].substitutions


def test_nonsmooth_abs_uses_jensen_graph_qe():
    x = sp.symbols("x", real=True)
    result = function_convexity(sp.Abs(x), (x,), return_result=True)
    assert result.classification == "convex"
    assert result.convex is True
    assert "epigraph_set_convexity" in result.method


def test_nonconvex_domain_is_reported_separately():
    x = sp.symbols("x", real=True)
    domain = sp.Or(x <= -1, x >= 1)
    result = function_convexity(x**2, (x,), domain=domain, return_result=True)
    assert result.classification == "nonconvex_domain"
    assert result.domain_convex is False


def test_natural_domain_is_intersected_automatically():
    x = sp.symbols("x", real=True)
    result = function_convexity(sp.sqrt(x), (x,), return_result=True)
    assert result.domain.has(x >= 0)
    assert result.classification == "strictly_concave"


def test_parameterized_global_quadratic_uses_hessian_stratification():
    x, a = sp.symbols("x a", real=True)
    result = function_convexity(a * x**2, (x,), parameters=(a,))
    assert isinstance(result, ParameterStratifiedResult)
    assert result.method == "function_convexity_hessian_parameter_stratification"
    assert result.select({a: 2}) == "strongly_convex"
    assert result.select({a: -2}) == "strongly_concave"
    assert result.select({a: 0}) == "affine"


def test_unsupported_transcendental_function_is_unknown_not_guessed():
    x = sp.symbols("x", real=True)
    result = function_convexity(sp.exp(x), (x,), return_result=True)
    assert result.classification == "unknown"
    assert result.convex is None
    assert result.concave is None
