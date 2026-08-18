import sympy as sp

from semialg.algebraic.rational_univariate import solve_formula_with_rur
from semialg.solve.find_instance import find_instance_text
from semialg.solve.reduce import reduce_text


def _truth(expr):
    value = sp.simplify(expr)
    if value == sp.true or value is sp.true:
        return True
    if value == sp.false or value is sp.false:
        return False
    return bool(value)


def test_rur_solves_and_filters_formula_branches():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x**2 + y**2, 1), sp.Eq(x - y, 0), x > 0, evaluate=False)

    result = solve_formula_with_rur(formula, (x, y))

    assert result is not None
    assert result.status == "satisfied"
    assert result.assignments == ({x: sp.sqrt(2) / 2, y: sp.sqrt(2) / 2},)


def test_find_instance_uses_rur_before_sampling_or_cad_for_finite_system():
    x, y = sp.symbols("x y", real=True)

    result = find_instance_text(
        "exists x, y. (x^2 + y^2 == 1) & (x - y == 0) & (x > 0)",
        symbols={"x": x, "y": y},
        variable_order=[x, y],
        domain="reals",
        strategy="auto",
        use_preprocess=False,
        return_result=True,
    )

    assert result.method == "rational_univariate_instance"
    assert result.found
    assert _truth(sp.And(sp.Eq(x**2 + y**2, 1), sp.Eq(x - y, 0), x > 0).subs(result.first()))


def test_reduce_uses_rur_and_projects_finite_existential_solutions():
    x, y = sp.symbols("x y", real=True)

    result = reduce_text(
        "exists x. (x^2 + y^2 == 1) & (x - y == 0)",
        symbols={"x": x, "y": y},
        variable_order=[x, y],
        strategy="auto",
        use_preprocess=False,
        return_result=True,
    )

    assert result.method == "rational_univariate"
    assert x not in result.result.free_symbols
    assert _truth(result.result.subs(y, sp.sqrt(2) / 2))
    assert _truth(result.result.subs(y, -sp.sqrt(2) / 2))
    assert not _truth(result.result.subs(y, 0))


def test_virtual_substitution_handles_pure_universal_quadratic_block():
    x, y = sp.symbols("x y", real=True)

    result = reduce_text(
        "forall x. x^2 + y^2 >= 0",
        symbols={"x": x, "y": y},
        variable_order=[x, y],
        strategy="auto",
        use_preprocess=False,
        return_result=True,
    )

    assert result.method == "quadratic_virtual_substitution"
    assert result.result == sp.true or sp.simplify(result.result) == sp.true
