import pytest
import sympy as sp

from semialg.function_graph import semialgebraic_formula_graph, semialgebraic_function_graph

x, target = sp.symbols("x target", real=True)


def _feasible(formula):
    # Eliminate unconditional affine equations before the exact linear oracle.
    # This also discards inactive branches containing nonlinear value graphs.
    equations = [
        atom.lhs - atom.rhs
        for atom in sp.And.make_args(formula)
        if isinstance(atom, sp.Equality) and sp.total_degree(atom.lhs - atom.rhs) <= 1
    ]
    if equations:
        variables = sorted(
            set().union(*(eq.free_symbols for eq in equations)), key=sp.default_sort_key
        )
        solutions = sp.linsolve(equations, variables)
        if solutions is sp.EmptySet:
            return False
        formula = formula.subs(
            dict(zip(variables, next(iter(solutions)), strict=True)), simultaneous=True
        )
    return sp.satisfiable(formula, use_lra_theory=True) is not False


@pytest.mark.parametrize("connective", [sp.Not, sp.Implies, sp.Equivalent, sp.Xor])
def test_boolean_graph_truth(connective):
    atom = sp.Abs(x) > 0
    formula = (
        connective(atom, evaluate=False)
        if connective is sp.Not
        else connective(atom, x < 0, evaluate=False)
    )
    graph = semialgebraic_formula_graph(formula)
    for value in (-1, 0, 1):
        expected = bool(formula.subs(x, value))
        assert _feasible(graph.formula.subs(x, value)) == expected


def test_piecewise_condition_complement():
    expression = sp.Piecewise((1, sp.And(sp.Abs(x) > 1, x > 0)), (2, True))
    graph = semialgebraic_function_graph(expression, target)
    for value in (0, 2):
        for output in (1, 2):
            condition = graph.formula.subs({x: value, target: output})
            assert _feasible(condition) == (expression.subs(x, value) == output)


@pytest.mark.parametrize("zero_value,point,expected", [(sp.sqrt(x), -1, 0), (1 / (x - 1), 1, 1)])
def test_heaviside_unused_zero_value(zero_value, point, expected):
    expression = sp.Heaviside(x, zero_value)
    assert expression.subs(x, point) == expected
    graph = semialgebraic_function_graph(expression, target)
    assert _feasible(graph.formula.subs({x: point, target: expected}))
    assert not _feasible(graph.formula.subs({x: point, target: expected + 1}))
