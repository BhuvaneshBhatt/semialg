import sympy as sp

from semialg.optimization import _graph_formula_for_expression
from semialg.optimization_active_sets import kkt_system


def test_kkt_multipliers_cannot_collide_with_user_symbol_names():
    x = sp.Symbol("x", real=True)
    user_lambda = sp.Symbol("semialg_lambda_0", real=True)
    equations, multipliers = kkt_system(x**2 + user_lambda * x, (x,), (x - 1,))

    assert len(multipliers) == 1
    assert isinstance(multipliers[0], sp.Dummy)
    assert multipliers[0] != user_lambda
    assert user_lambda in equations[-1].free_symbols


def test_graph_auxiliaries_are_collision_free_dummies():
    x = sp.Symbol("x", real=True)
    target = sp.Symbol("y", real=True)
    user_aux = sp.Symbol("semialg_graph_aux", real=True)

    formula, auxiliaries = _graph_formula_for_expression(sp.Abs(x) + user_aux, target)

    assert auxiliaries
    assert all(isinstance(aux, sp.Dummy) for aux in auxiliaries)
    assert user_aux in formula.free_symbols
    assert all(aux != user_aux for aux in auxiliaries)
