import sympy as sp

from semialg import integrate_over_region, semialgebraic_measure
from semialg.conditional import ParameterStratifiedResult


def test_parametric_interval_integral_covers_empty_and_feasible_fibers():
    x, a = sp.symbols("x a", real=True)

    result = integrate_over_region(
        x,
        sp.And(x >= 0, x <= a),
        [x],
        parameters=[a],
        return_stratified=True,
    )

    assert isinstance(result, ParameterStratifiedResult)
    assert result.complete is True
    assert result.coverage_condition is sp.true
    assert sp.simplify(result.select({a: 2}) - 2) == 0
    assert result.select({a: -1}) == 0
    assert result.select({a: 0}) == 0


def test_parametric_measure_returns_exact_piecewise_length():
    x, a = sp.symbols("x a", real=True)

    result = semialgebraic_measure(
        sp.And(x >= 0, x <= a),
        [x],
        parameters=[a],
        return_stratified=True,
    )

    assert isinstance(result, ParameterStratifiedResult)
    assert result.select({a: 3}) == 3
    assert result.select({a: -2}) == 0


def test_parametric_box_integral_in_two_dimensions():
    x, y, a = sp.symbols("x y a", real=True)

    result = integrate_over_region(
        x + y,
        sp.And(x >= 0, x <= a, y >= 0, y <= 1),
        [x, y],
        parameters=[a],
        return_stratified=True,
    )

    assert sp.simplify(result.select({a: 2}) - 3) == 0
    assert result.select({a: -1}) == 0


def test_parameters_require_stratified_result_for_region_integrals():
    x, a = sp.symbols("x a", real=True)

    try:
        integrate_over_region(1, sp.And(x >= 0, x <= a), [x], parameters=[a])
    except ValueError as exc:
        assert "return_stratified=True" in str(exc)
    else:
        raise AssertionError("expected ValueError")
