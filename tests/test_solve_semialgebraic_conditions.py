import pytest
import sympy as sp

from semialg import solve_semialgebraic
from semialg.parameters import solvability_conditions


def test_output_conditions_for_parameterized_quadratic_equality():
    x, a, b = sp.symbols("x a b", real=True)

    condition = solve_semialgebraic(
        [sp.Eq(x**2 + a * x + b, 0)],
        [x],
        parameters=[a, b],
        count=0,
        output="conditions",
    )

    assert condition == (a**2 - 4 * b >= 0)


def test_output_conditions_for_parameterized_quadratic_inequality():
    x, a = sp.symbols("x a", real=True)

    condition = solve_semialgebraic(
        [x**2 + a < 0],
        [x],
        parameters=[a],
        count=0,
        output="conditions",
    )

    assert condition == (a < 0)


def test_output_conditions_include_degenerate_linear_equation_case():
    x, a, b = sp.symbols("x a b", real=True)

    unconditional = solve_semialgebraic(
        sp.Eq(a * x, 0), [x], parameters=[a], count=0, output="conditions"
    )
    conditional = solve_semialgebraic(
        sp.Eq(a * x + b, 0), [x], parameters=[a, b], count=0, output="conditions"
    )

    assert unconditional is sp.true
    assert conditional.subs({a: 0, b: 0}) is sp.true
    assert conditional.subs({a: 0, b: 1}) is sp.false
    assert conditional.subs({a: 1, b: 1}) is sp.true


@pytest.mark.parametrize("solver", [solve_semialgebraic, solvability_conditions])
def test_explicit_parameter_list_rejects_undeclared_symbols(solver):
    x, a, b = sp.symbols("x a b", real=True)

    with pytest.raises(ValueError, match=r"undeclared symbolic parameters: b"):
        solver(sp.Eq(b * x + 1, 0), [x], parameters=[a])


def test_string_variables_and_parameters_preserve_formula_symbol_identity():
    x = sp.Symbol("x", real=True)
    a = sp.Symbol("a")

    result = solve_semialgebraic(sp.Eq(x, a), ["x"], parameters=["a"], count=0)

    assert result.variables == (x,)
    assert result.parameters == (a,)


def test_result_object_still_records_parameter_conditions():
    x, a, b = sp.symbols("x a b", real=True)

    result = solve_semialgebraic(
        [sp.Eq(x**2 + a * x + b, 0)],
        [x],
        parameters=[a, b],
        count=0,
    )

    assert result.parameter_conditions == (a**2 - 4 * b >= 0)
    assert result.parameters == (a, b)
    assert result.diagnostics["capabilities"]["parameter_conditions"] is True


def test_output_selectors_reuse_solution_result():
    x = sp.symbols("x", real=True)

    components = solve_semialgebraic([x**2 <= 1], [x], count=0, output="components")
    formula = solve_semialgebraic([x**2 <= 1], [x], count=0, output="formula")
    satisfiable = solve_semialgebraic([x**2 <= 1], [x], count=0, output="satisfiable")

    assert len(components) == 1
    assert formula == sp.And(x >= -1, x <= 1)
    assert satisfiable is True
