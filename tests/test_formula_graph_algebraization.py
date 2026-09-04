import sympy as sp

from semialg import function_range
from semialg.function_graph import semialgebraic_formula_graph


def test_formula_graph_algebraizes_abs_in_constraint() -> None:
    x, a = sp.symbols("x a", real=True)
    graph = semialgebraic_formula_graph(sp.Abs(x) <= a)

    assert not graph.formula.has(sp.Abs)
    assert graph.auxiliary_variables


def test_function_range_accepts_semialgebraic_function_in_constraint() -> None:
    x, t = sp.symbols("x t", real=True)
    result = function_range(x, sp.Abs(x) <= 2, [x], value_symbol=t)

    assert sp.simplify(result.subs(t, 0)) is sp.true
    assert sp.simplify(result.subs(t, 2)) is sp.true
    assert sp.simplify(result.subs(t, 3)) is sp.false


def test_piecewise_branch_condition_can_contain_semialgebraic_function() -> None:
    x, t = sp.symbols("x t", real=True)
    expr = sp.Piecewise((x, sp.Abs(x) <= 1), (sp.Integer(2), True))
    result = function_range(expr, sp.And(x >= -2, x <= 2), [x], value_symbol=t)

    for value in (-1, 0, 1, 2):
        assert sp.simplify(result.subs(t, value)) is sp.true
    assert sp.simplify(result.subs(t, sp.Rational(3, 2))) is sp.false
