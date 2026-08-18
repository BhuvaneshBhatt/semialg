import sympy as sp

from semialg import function_range


def test_function_range_closed_interval_metadata():
    x, t = sp.symbols("x t", real=True)
    result = function_range(x**2, sp.And(x >= -2, x <= 3), [x], value_symbol=t, return_result=True)

    assert sp.simplify(result.range_condition ^ sp.And(t >= 0, t <= 9)) is sp.false
    assert result.lower_bound == 0
    assert result.upper_bound == 9
    assert result.lower_bound_attained is True
    assert result.upper_bound_attained is True
    assert result.is_interval is True
    assert result.interval_count == 1


def test_function_range_open_interval_metadata():
    x, t = sp.symbols("x t", real=True)
    result = function_range(
        2 * x + 1, sp.And(x > 0, x < 1), [x], value_symbol=t, return_result=True
    )

    assert sp.simplify(result.range_condition ^ sp.And(t > 1, t < 3)) is sp.false
    assert result.lower_bound == 1
    assert result.upper_bound == 3
    assert result.lower_bound_attained is False
    assert result.upper_bound_attained is False
    assert result.is_interval is True
    assert result.interval_count == 1


def test_function_range_half_unbounded_metadata():
    x, t = sp.symbols("x t", real=True)
    result = function_range(1 / x, x > 0, [x], value_symbol=t, return_result=True)

    assert sp.simplify(result.range_condition ^ (t > 0)) is sp.false
    assert result.lower_bound == 0
    assert result.upper_bound == sp.oo
    assert result.lower_bound_attained is False
    assert result.upper_bound_attained is False
    assert result.is_interval is True
    assert result.interval_count == 1


def test_function_range_disconnected_metadata():
    x, t = sp.symbols("x t", real=True)
    result = function_range(x, sp.Or(x <= -1, x >= 1), [x], value_symbol=t, return_result=True)

    assert sp.simplify(result.range_condition ^ (sp.Abs(t) >= 1)) is sp.false
    assert result.lower_bound == -sp.oo
    assert result.upper_bound == sp.oo
    assert result.lower_bound_attained is False
    assert result.upper_bound_attained is False
    assert result.is_interval is False
    assert result.interval_count == 2
