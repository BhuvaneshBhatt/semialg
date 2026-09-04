import pytest
import sympy as sp

from semialg import is_equal
from semialg.algebraic.rational_univariate import solve_zero_dimensional_system_with_rur
from semialg.formula import parse_quant_form_text
from semialg.qe.complete import qe_by_complete_cad
from semialg.qe.virtual_substitution import try_quadratic_virtual_substitution_qe
from semialg.solve import find_instance

x, y, a = sp.symbols("x y a", real=True)


@pytest.mark.parametrize(
    "text",
    [
        "exists x. x^2 <= 1",
        "exists x. x^2 + a <= 1",
        "exists x. (x^2 <= a) & (x >= 0)",
        "exists x. (x^2 == a) & (x >= 0)",
        "exists x. (x^2 + x <= a)",
    ],
)
def test_virtual_substitution_matches_complete_cad(text):
    parsed = parse_quant_form_text(text, symbols={"x": x, "a": a}, variable_order=[a, x])
    vs = try_quadratic_virtual_substitution_qe(parsed.vars, parsed.quantifiers, parsed.matrix_expr)
    assert vs is not None and vs.status == "complete"
    cad = qe_by_complete_cad(
        parsed.vars,
        parsed.quantifiers,
        parsed.matrix,
        free_variables=[a] if a in parsed.matrix_expr.free_symbols else [],
        variable_order_strategy="preserve",
        use_presolve=False,
        return_result=True,
    )
    free = [a] if a in vs.formula.free_symbols or a in cad.formula.free_symbols else []
    if free:
        assert is_equal(vs.formula, cad.formula, free)
    else:
        assert bool(vs.formula) == bool(cad.formula)


@pytest.mark.parametrize(
    "equations,expected_count",
    [
        ([x**2 - 1], 2),
        ([x**2 - 4], 2),
        ([x**3 - x], 3),
        ([x**2 + 1], 0),
        ([x**2 - 1, y - x], 2),
        ([x**2 + y**2 - 1, x - y], 2),
    ],
)
def test_rur_real_solution_count_matches_cad_existence(equations, expected_count):
    variables = [x] if all(y not in eq.free_symbols for eq in equations) else [x, y]
    roots = solve_zero_dimensional_system_with_rur(equations, variables, real=True)
    assert len(roots) == expected_count
    formula = sp.And(*(sp.Eq(eq, 0) for eq in equations))
    witness = find_instance(formula, variables, return_result=False)
    assert (witness is not None) == (expected_count > 0)


@pytest.mark.parametrize(
    "formula,variables",
    [
        (sp.And(x >= 0, x <= 1), [x]),
        (sp.And(x**2 <= 4, x > 1), [x]),
        (sp.And(x**2 + y**2 < 1, x > 0), [x, y]),
        (sp.And(x >= -1, x <= 1, y >= -1, y <= 1), [x, y]),
        (sp.And(x > 1, x < 0), [x]),
    ],
)
def test_find_instance_agrees_with_existential_complete_qe(formula, variables):
    witness = find_instance(formula, variables, return_result=False)
    from semialg.formula import parse_formula

    result = qe_by_complete_cad(
        variables,
        [("exists", v) for v in variables],
        parse_formula(formula),
        free_variables=[],
        variable_order_strategy="preserve",
        return_result=True,
    )
    assert (witness is not None) == bool(result.truth_value)
