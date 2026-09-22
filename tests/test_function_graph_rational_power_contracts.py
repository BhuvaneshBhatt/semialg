import sympy as sp

from semialg import is_satisfiable
from semialg.function_graph import semialgebraic_function_graph


def test_negative_rational_power_excludes_zero_and_is_exact():
    x, y = sp.symbols("x y", real=True)
    graph = semialgebraic_function_graph(x ** sp.Rational(-2, 3), y).formula
    assert is_satisfiable(sp.And(graph, sp.Eq(x, 8), sp.Eq(y, sp.Rational(1, 4))), [x, y])
    assert not is_satisfiable(sp.And(graph, sp.Eq(x, 0)), [x, y])


def test_explicit_odd_real_root_accepts_negative_base():
    x, y = sp.symbols("x y", real=True)
    expr = sp.real_root(x, 3)
    graph = semialgebraic_function_graph(expr, y).formula
    assert is_satisfiable(sp.And(graph, sp.Eq(x, -8), sp.Eq(y, -2)), [x, y])
