import sympy as sp

from semialg import equivalent, simplify_boole, simplify_piecewise


def same_logic(lhs, rhs, variables):
    return equivalent(lhs, rhs, variables)


def test_simplify_boole_removes_redundant_conjunct():
    x = sp.Symbol("x", real=True)
    simplified = simplify_boole(sp.And(x > 0, x >= 0), [x])
    assert same_logic(simplified, x > 0, [x])


def test_simplify_boole_detects_tautology_and_unsat():
    x, y = sp.symbols("x y", real=True)
    assert simplify_boole(sp.Or(x < 0, x >= 0), [x]) == sp.true
    assert simplify_boole(sp.And(x**2 < 0, y > 0), [x, y]) == sp.false


def test_simplify_boole_removes_redundant_disjunct():
    x, y = sp.symbols("x y", real=True)
    expr = sp.Or(x > 0, sp.And(x > 1, y > 0))
    simplified = simplify_boole(expr, [x, y])
    assert same_logic(simplified, x > 0, [x, y])


def test_simplify_piecewise_removes_impossible_branch():
    x = sp.Symbol("x", real=True)
    expr = sp.Piecewise((1, x**2 < 0), (2, x > 0), (3, True))
    simplified = simplify_piecewise(expr, [x])
    expected = sp.Piecewise((2, x > 0), (3, True), evaluate=False)
    assert simplified == expected


def test_simplify_piecewise_collapses_under_assumptions():
    x = sp.Symbol("x", real=True)
    expr = sp.Piecewise((x, x >= 0), (-x, True))
    assert simplify_piecewise(expr, [x], assumptions=x >= 0) == x
    assert simplify_piecewise(expr, [x], assumptions=x <= 0) == -x


def test_simplify_piecewise_merges_equal_adjacent_values():
    x = sp.Symbol("x", real=True)
    expr = sp.Piecewise((1, x < 0), (1, x >= 0), (2, True), evaluate=False)
    assert simplify_piecewise(expr, [x]) == 1
