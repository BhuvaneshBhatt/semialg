import sympy as sp

from semialg import is_satisfiable, quantifier_eliminate


def test_boolean_contradiction_avoids_algebraic_solver():
    x = sp.symbols("x", real=True)
    relation = x**7 - x + 1 > 0
    other = x**5 + 2 > 0
    formula = sp.And(relation, sp.Or(sp.Not(relation), other), sp.Not(other), evaluate=False)
    result = is_satisfiable(formula, [x], return_result=True)
    assert result.satisfiable is False
    assert result.method == "propositional"


def test_vacuous_quantifier_is_removed():
    x, y = sp.symbols("x y", real=True)
    result = quantifier_eliminate(x > 0, [("exists", y)], return_result=True)
    assert result.formula == (x > 0)
    assert result.method == "quantifier-free"
    assert result.quantified_variables == ()


def test_large_boolean_skeleton_does_not_force_exponential_minimization():
    xs = sp.symbols("x0:10", real=True)
    atoms = tuple(x > 0 for x in xs)
    formula = sp.Or(*atoms, evaluate=False)
    result = is_satisfiable(formula, xs, return_result=True)
    assert result.satisfiable is True
    assert result.method != "propositional"
