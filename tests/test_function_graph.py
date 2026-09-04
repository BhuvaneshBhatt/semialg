import sympy as sp

from semialg.function_graph import semialgebraic_function_graph


def test_odd_real_root_integer_powers_use_direct_exact_graph():
    x, target = sp.symbols("x target", real=True)
    for integer_power in (1, 2, 4, -1):
        expression = sp.real_root(x, 3) ** integer_power
        graph = semialgebraic_function_graph(expression, target)
        expected_diagnostics = (
            ("explicit_real_root_semantics",)
            if integer_power == 1
            else ("explicit_real_root_power_semantics",)
        )
        assert graph.diagnostics == expected_diagnostics
        assert len(graph.auxiliary_variables) == 1
        base_auxiliary = graph.auxiliary_variables[0]
        for base_value in (-8, -1, 1, 8):
            expected = sp.simplify(expression.subs(x, base_value))
            assignment = {x: base_value, target: expected, base_auxiliary: base_value}
            assert sp.simplify(graph.formula.subs(assignment)) is sp.true


def test_heaviside_has_exact_semialgebraic_graph_including_zero_value():
    x, target = sp.symbols("x target", real=True)
    graph = semialgebraic_function_graph(sp.Heaviside(x, sp.Rational(1, 3)), target)

    for x_value, expected in ((-2, 0), (0, sp.Rational(1, 3)), (5, 1)):
        reduced = sp.simplify(graph.formula.subs({x: x_value, target: expected}))
        from semialg import is_satisfiable

        assert is_satisfiable(reduced, graph.auxiliary_variables)
