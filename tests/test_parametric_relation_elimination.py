import sympy as sp

from semialg import function_range, semialgebraic_minimize
from semialg.optimization_results import ParametricFunctionRangeResult, ParametricOptimizationResult


def test_function_range_can_optionally_eliminate_parametric_relation():
    x, a, t = sp.symbols("x a t", real=True)

    result = function_range(
        a,
        sp.Eq(x, 0),
        [x],
        value_symbol=t,
        parameters=[a],
        return_stratified=True,
        eliminate_quantifiers=True,
    )

    value = result.branches[0].value
    assert isinstance(value, ParametricFunctionRangeResult)
    assert value.quantifier_free is True
    assert value.quantifiers == ()
    assert not any(
        symbol.name.startswith("_semialg_range_") for symbol in value.formula.free_symbols
    )


def test_optimization_can_optionally_eliminate_parametric_relation():
    x, a = sp.symbols("x a", real=True)

    result = semialgebraic_minimize(
        a,
        sp.Eq(x, 0),
        [x],
        parameters=[a],
        return_stratified=True,
        eliminate_quantifiers=True,
        return_result=True,
    )

    value = result.branches[0].value
    assert isinstance(value, ParametricOptimizationResult)
    assert value.quantifier_free is True
    assert value.quantifiers == ()


def test_relation_elimination_is_opt_in():
    x, a, t = sp.symbols("x a t", real=True)

    result = function_range(
        a,
        sp.Eq(x, 0),
        [x],
        value_symbol=t,
        parameters=[a],
        return_stratified=True,
    )

    value = result.branches[0].value
    assert value.quantifier_free is False
    assert value.quantifiers


def test_eliminate_quantifiers_requires_stratified_mode():
    x = sp.symbols("x", real=True)

    try:
        function_range(x, x >= 0, [x], eliminate_quantifiers=True)
    except ValueError as exc:
        assert "return_stratified=True" in str(exc)
    else:
        raise AssertionError("expected ValueError")
