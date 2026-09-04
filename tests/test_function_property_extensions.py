import sympy as sp

from semialg import (
    function_convexity,
    function_mapping_properties,
    function_smoothness,
)


def test_convexity_strong_and_strict_primary_results():
    x = sp.symbols("x", real=True)
    assert function_convexity(x**2, [x]) == "strongly_convex"
    assert function_convexity(x**4, [x]) == "strictly_convex"
    assert function_convexity(-(x**2), [x]) == "strongly_concave"


def test_parameterized_quadratic_strong_curvature():
    x, a = sp.symbols("x a", real=True)
    result = function_convexity(a * x**2, [x])
    assert result.select({a: 2}) == "strongly_convex"
    assert result.select({a: 0}) == "affine"
    assert result.select({a: -2}) == "strongly_concave"


def test_all_convexity_properties_and_pseudoconvex_stationary_logic():
    x = sp.symbols("x", real=True)
    quartic = function_convexity(x**4, [x], properties="all", return_result=True)
    assert quartic.strictly_convex is True
    assert quartic.quasiconvex is True
    assert quartic.strictly_quasiconvex is True
    assert quartic.pseudoconvex is True

    cubic = function_convexity(x**3, [x], properties="all", return_result=True)
    assert cubic.quasiconvex is True
    assert cubic.pseudoconvex is False

    absolute = function_convexity(sp.Abs(x), [x], properties="all", return_result=True)
    assert absolute.convex is True
    assert absolute.pseudoconvex is False


def test_algebraic_log_convexity_and_concavity():
    x = sp.symbols("x", real=True)
    result = function_convexity(
        1 / (x + 1), [x], domain=x > -1, properties="all", return_result=True
    )
    assert result.log_convex is True

    result = function_convexity(
        1 - x**2,
        [x],
        domain=sp.And(x > -1, x < 1),
        properties="all",
        return_result=True,
    )
    assert result.log_concave is True


def test_function_smoothness_polynomial_abs_sign_piecewise():
    x = sp.symbols("x", real=True)
    polynomial = function_smoothness(x**3, x)
    assert polynomial.continuous is True
    assert polynomial.smooth is True
    assert polynomial.differentiability_order == sp.oo

    absolute = function_smoothness(sp.Abs(x), x)
    assert absolute.continuous is True
    assert absolute.differentiability_order == 0
    assert absolute.smooth is False
    assert absolute.derivative_exceptions[1] == sp.Eq(x, 0)

    sign = function_smoothness(sp.sign(x), x)
    assert sign.continuous is False
    assert sign.continuity_exceptions == sp.Eq(x, 0)

    piecewise = sp.Piecewise((x**2, x < 0), (2 * x**2, True))
    joined = function_smoothness(piecewise, x)
    assert joined.continuous is True
    assert joined.differentiability_order == 1
    assert joined.derivative_exceptions[2] == sp.Eq(x, 0)


def test_function_mapping_properties_scalar_maps():
    x = sp.symbols("x", real=True)
    identity = function_mapping_properties(x, x)
    assert identity.injective is True
    assert identity.surjective is True
    assert identity.bijective is True

    square = function_mapping_properties(x**2, x)
    assert square.injective is False
    assert square.surjective is False
    assert square.bijective is False
    assert square.collision_witness is not None
    assert square.missing_value_witness is not None

    y0 = sp.Symbol("y0", real=True)
    nonnegative = function_mapping_properties(x**2, x, codomain=y0 >= 0)
    assert nonnegative.injective is False
    assert nonnegative.surjective is True


def test_function_mapping_properties_vector_affine_map():
    x, y = sp.symbols("x y", real=True)
    result = function_mapping_properties((x + y, x - y), (x, y))
    assert result.injective is True
    assert result.surjective is True
    assert result.bijective is True
