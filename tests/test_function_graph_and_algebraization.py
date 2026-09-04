import sympy as sp

from semialg import equivalent, function_domain, function_range
from semialg._function_algebraization import (
    algebraic_back_substitute,
    exact_algebraize_function_problem,
)
from semialg.function_graph import semialgebraic_formula_graph, semialgebraic_function_graph


def test_principal_pow_and_explicit_real_root_have_distinct_real_domains() -> None:
    x = sp.symbols("x", real=True)

    assert function_domain(x ** sp.Rational(1, 3), [x]) == (x >= 0)
    assert function_domain(x ** sp.Rational(-1, 3), [x]) == (x > 0)
    assert function_domain(sp.real_root(x, 3), [x]) is sp.true


def test_shared_graph_preserves_real_root_semantics() -> None:
    x, y = sp.symbols("x y", real=True)
    graph = semialgebraic_function_graph(sp.real_root(x, 3), y)

    # The canonical real-root graph is exactly y**3 == x after eliminating
    # internal auxiliaries; range/domain projections exercise that fact.
    assert graph.diagnostics == ("explicit_real_root_semantics",)
    assert function_range(sp.real_root(x, 3), variables=[x], value_symbol=y) is sp.true


def test_nested_principal_radicals_use_graph_projection_for_domain() -> None:
    x = sp.symbols("x", real=True)
    result = function_domain(sp.sqrt(1 - sp.sqrt(x)), [x])
    assert equivalent(result, sp.And(x >= 0, x <= 1), [x])


def test_algebraic_back_substitution_removes_uniquely_solved_variable() -> None:
    x, y = sp.symbols("x y", real=True)
    result = algebraic_back_substitute(sp.exp(y), sp.Eq(y, 2 * x), (x, y))

    assert result.expression == sp.exp(2 * x)
    assert result.constraints is sp.true
    assert result.variables == (x,)
    assert result.substitutions == ((y, 2 * x),)


def test_commensurate_trigonometric_range_is_algebraized_exactly() -> None:
    x, t = sp.symbols("x t", real=True)
    result = function_range(
        sp.sin(x) + sp.cos(2 * x), variables=[x], value_symbol=t, return_result=True
    )

    assert result.method == "algebraized_commensurate_trigonometric"
    assert equivalent(result.range_condition, sp.And(t >= -2, t <= sp.Rational(9, 8)), [t])


def test_trigonometric_constraint_is_algebraized_with_same_phase_variable() -> None:
    x, t = sp.symbols("x t", real=True)
    result = function_range(sp.sin(x), sp.cos(x) >= 0, [x], value_symbol=t, return_result=True)

    assert result.method == "algebraized_commensurate_trigonometric"
    assert equivalent(result.range_condition, sp.And(t >= -1, t <= 1), [t])


def test_commensurate_exponential_and_hyperbolic_ranges() -> None:
    x, t = sp.symbols("x t", real=True)

    exp_result = function_range(
        sp.exp(x) + sp.exp(-x), variables=[x], value_symbol=t, return_result=True
    )
    cosh_result = function_range(sp.cosh(2 * x), variables=[x], value_symbol=t, return_result=True)

    assert exp_result.method == "algebraized_commensurate_exponential_hyperbolic"
    assert cosh_result.method == "algebraized_commensurate_exponential_hyperbolic"
    assert equivalent(exp_result.range_condition, t >= 2, [t])
    assert equivalent(cosh_result.range_condition, t >= 1, [t])


def test_back_substitution_precedes_exponential_algebraization() -> None:
    x, y, t = sp.symbols("x y t", real=True)
    result = function_range(
        sp.exp(y) + sp.exp(-y),
        sp.Eq(y, 2 * x),
        [x, y],
        value_symbol=t,
        return_result=True,
    )

    assert result.method == "algebraized_commensurate_exponential_hyperbolic"
    assert equivalent(result.range_condition, t >= 2, [t])
    assert result.diagnostics["back_substitutions"] == (("y", "2*x"),)


def test_noncommensurate_trig_frequencies_are_not_relaxed_to_independent_circles() -> None:
    x = sp.symbols("x", real=True)
    algebraized = exact_algebraize_function_problem(
        sp.sin(x) + sp.sin(sp.sqrt(2) * x), sp.true, (x,)
    )
    assert algebraized is None


def test_trig_algebraization_declines_when_bare_source_variable_remains() -> None:
    x = sp.symbols("x", real=True)
    algebraized = exact_algebraize_function_problem(sp.sin(x) + x, sp.true, (x,))
    assert algebraized is None


def test_formula_graph_supports_xor_implies_and_equivalent():
    x = sp.Symbol("x", real=True)
    formula = sp.Equivalent(sp.Implies(sp.Abs(x) <= 1, x <= 2), sp.Xor(x < 0, x >= 0))
    graph = semialgebraic_formula_graph(formula)
    assert graph.formula is not None
    assert graph.auxiliary_variables
