import sympy as sp

from semialg.algebraic.rational_univariate import solve_formula_with_rur
from semialg.solve import find_instance_text, reduce_text


def test_rur_distributes_small_boolean_branches_inside_conjunction():
    x, y = sp.symbols("x y")
    result = solve_formula_with_rur(
        sp.And(sp.Eq(y, 2), sp.Or(sp.Eq(x, 0), sp.Eq(x, 1))),
        (x, y),
    )

    assert result is not None
    assert result.points == ((0, 2), (1, 2))


def test_rur_ignores_tautological_equalities_and_detects_false_numeric_equalities():
    x = sp.symbols("x")
    sat = solve_formula_with_rur(sp.And(sp.Eq(x, 1), sp.Eq(0, 0)), (x,))
    unsat = solve_formula_with_rur(sp.And(sp.Eq(x, 1), sp.Eq(0, 1)), (x,))

    assert sat is not None
    assert sat.points == ((1,),)
    assert unsat is not None
    assert unsat.status == "unsat"


def test_reduce_rur_projects_finite_existential_solution_to_free_variable():
    x, y = sp.symbols("x y", real=True)

    result = reduce_text(
        "exists x. (x^2 + y^2 == 1) & (x - y == 0)",
        symbols={"x": x, "y": y},
        variable_order=[y, x],
        strategy="auto",
        use_preprocess=False,
        return_result=True,
    )

    assert result.method == "rational_univariate"
    assert x not in result.result.free_symbols
    assert sp.simplify(result.result.subs(y, sp.sqrt(2) / 2)) == sp.true
    assert sp.simplify(result.result.subs(y, -sp.sqrt(2) / 2)) == sp.true
    assert sp.simplify(result.result.subs(y, 0)) == sp.false


def test_virtual_substitution_reconstructs_witness_after_rur_edge_case():
    x, y = sp.symbols("x y", real=True)

    result = find_instance_text(
        "exists x. (x > 0) & (x^2 + y^2 <= 1) & (y == 1/2)",
        symbols={"x": x, "y": y},
        variable_order=[x, y],
        strategy="auto",
        use_preprocess=False,
        return_result=True,
    )

    assert result.method == "quadratic_virtual_substitution_instance"
    witness = result.first()
    assert witness is not None
    assert sp.simplify(witness[y] - sp.Rational(1, 2)) == 0
    assert sp.simplify((witness[x] > 0) & (witness[x] ** 2 + witness[y] ** 2 <= 1)) == sp.true
