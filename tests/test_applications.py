import pytest
import sympy as sp

from semialg.applications import (
    exact_optimization_benchmark,
    robust_parameter_analysis,
    robust_parameter_region,
    validate_formula_equivalence,
    validate_identity,
    validate_numeric_optimization,
    validate_range,
)


def test_robust_parameter_region_for_quadratic_nonnegativity():
    x, a, b = sp.symbols("x a b", real=True)
    condition = robust_parameter_region(x**2 + a * x + b >= 0, [x], [a, b])
    assert sp.simplify_logic(sp.Equivalent(condition, a**2 - 4 * b <= 0)) is sp.true


def test_robust_parameter_analysis_distinguishes_feasible_and_robust():
    x, a = sp.symbols("x a", real=True)
    result = robust_parameter_analysis(
        x**2 <= a,
        [x],
        [a],
        operating_domain=sp.And(x >= -1, x <= 1),
    )
    assert sp.simplify_logic(sp.Equivalent(result.feasible_condition, a >= 0)) is sp.true
    assert sp.simplify_logic(sp.Equivalent(result.robust_condition, a >= 1)) is sp.true


def test_robust_parameter_region_rejects_unknown_quantifier():
    x, a = sp.symbols("x a", real=True)
    with pytest.raises(ValueError, match="quantifier"):
        robust_parameter_region(x <= a, [x], [a], quantifier="sometimes")


def test_validate_identity_returns_counterexample_for_false_claim():
    x = sp.symbols("x", real=True)
    result = validate_identity(x**2, x, [x])
    assert not result.valid
    assert result.counterexample is not None


def test_validate_identity_under_assumptions():
    x = sp.symbols("x", real=True)
    result = validate_identity(sp.sqrt(x**2), x, [x], assumptions=x >= 0)
    assert result.valid


def test_validate_formula_equivalence():
    x = sp.symbols("x", real=True)
    result = validate_formula_equivalence(sp.And(x >= 0, x <= 1), sp.And(x <= 1, x >= 0), [x])
    assert result.valid


def test_validate_range():
    x, t = sp.symbols("x t", real=True)
    result = validate_range(
        x**2, sp.And(t >= 0, t <= 1), [x], constraints=sp.And(x >= -1, x <= 1), value_symbol=t
    )
    assert result.valid


def test_exact_optimization_benchmark_and_numeric_check():
    x = sp.symbols("x", real=True)
    benchmark = exact_optimization_benchmark(x**2, [x >= 1], [x])
    assert benchmark.certified
    assert benchmark.exact_value == 1
    check = validate_numeric_optimization(benchmark, 1.000001, atol=1e-5)
    assert check.within_tolerance
    assert check.absolute_error < 1e-5


def test_numeric_check_rejects_negative_tolerance():
    x = sp.symbols("x", real=True)
    benchmark = exact_optimization_benchmark(x**2, [x >= 1], [x])
    with pytest.raises(ValueError, match="nonnegative"):
        validate_numeric_optimization(benchmark, 1.0, atol=-1)


def test_robust_parameter_analysis_preserves_string_symbol_identity_and_rejects_overlap():
    x = sp.Symbol("x")
    a = sp.Symbol("a")
    result = robust_parameter_analysis(x <= a, ["x"], ["a"], operating_domain=x >= 0)
    assert result.variables == (x,)
    assert result.parameters == (a,)
    with pytest.raises(ValueError, match="disjoint"):
        robust_parameter_analysis(x >= 0, [x], [x])


def test_validation_rejects_false_formula_equivalence_and_false_range():
    x, t = sp.symbols("x t", real=True)
    formula = validate_formula_equivalence(x >= 0, x > 0, [x])
    assert not formula.valid
    assert formula.counterexample is not None
    range_check = validate_range(
        x**2,
        sp.And(t >= 0, t < 1),
        [x],
        constraints=sp.And(x >= -1, x <= 1),
        value_symbol=t,
    )
    assert not range_check.valid


def test_optimization_benchmark_supports_maximum_and_detects_bad_numeric_value():
    x = sp.Symbol("x", real=True)
    benchmark = exact_optimization_benchmark(x, [x >= 0, x <= 2], [x], kind="max")
    assert benchmark.certified
    assert benchmark.exact_value == 2
    check = validate_numeric_optimization(benchmark, 1.8, atol=1e-3)
    assert not check.within_tolerance


def test_optimization_benchmark_rejects_unknown_kind():
    x = sp.Symbol("x", real=True)
    with pytest.raises(ValueError, match="kind"):
        exact_optimization_benchmark(x**2, None, [x], kind="median")
