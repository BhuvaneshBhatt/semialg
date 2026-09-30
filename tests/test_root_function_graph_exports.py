import sympy as sp

import semialg


def test_semialgebraic_function_graph_is_supported_from_root_api():
    x, y = sp.symbols("x y", real=True)
    graph = semialg.semialgebraic_function_graph(sp.Abs(x), y)

    assert isinstance(graph, semialg.SemialgebraicFunctionGraph)
    assert graph.expression == sp.Abs(x)
    assert graph.target == y
    assert bool(graph.formula.subs({x: -2, y: 2}))


def test_unsupported_function_graph_is_supported_root_exception():
    x, y = sp.symbols("x y", real=True)
    try:
        semialg.semialgebraic_function_graph(sp.sin(x), y)
    except semialg.UnsupportedFunctionGraph:
        pass
    else:
        raise AssertionError("unsupported transcendental graph should raise")
