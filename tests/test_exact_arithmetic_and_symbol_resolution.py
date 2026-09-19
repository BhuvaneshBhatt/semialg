from __future__ import annotations

import pytest
import sympy as sp

import semialg
from semialg.domain_solve import normalize_domain_sensitive_constraints
from semialg.exact_arithmetic import compare_exact_reals, exact_truth
from semialg.reasoning import region_subset
from semialg.region_integrate import _finite_real_roots as integral_real_roots


def test_measure_string_variable_reuses_unassumed_formula_symbol():
    x = sp.Symbol("x")
    result = semialg.semialgebraic_measure(x > 0, ["x"], bounds={"x": (0, 1)}, return_result=True)
    assert result.value == 1
    assert result.variables == (x,)


def test_region_integral_string_variable_reuses_unassumed_formula_symbol():
    x = sp.Symbol("x")
    result = semialg.integrate_over_region(
        1, x > 0, ["x"], bounds={"x": (0, 1)}, return_result=True
    )
    assert result.value == 1
    assert result.variables == (x,)


def test_region_moment_string_variable_reuses_unassumed_formula_symbol():
    x = sp.Symbol("x")
    result = semialg.region_moment(
        x >= 0, ["x"], powers=(1,), bounds={"x": (0, 1)}, return_result=True
    )
    assert result.value == sp.Rational(1, 2)
    assert result.variables == (x,)


def test_domain_normalization_string_variable_reuses_unassumed_symbol():
    x = sp.Symbol("x")
    result = normalize_domain_sensitive_constraints(sp.sqrt(x) >= 0, ["x"])
    assert result.variables == (x,)
    assert len(result.variables) == 1


def test_reasoning_string_variable_does_not_duplicate_assumption_variant():
    x = sp.Symbol("x")
    assert region_subset(x > 1, x > 0, ["x"])


def test_symbolic_simplification_string_variable_preserves_input_symbol():
    x = sp.Symbol("x")
    result = semialg.simplify_boole(sp.And(x > 0, x >= 0), ["x"], return_result=True)
    assert result.variables == (x,)
    assert result.formula == (x > 0)


def test_exact_comparison_handles_close_algebraic_values_without_numeric_guess():
    lower = sp.Rational(1414213562373095, 10**15)
    assert compare_exact_reals(sp.sqrt(2), lower) > 0
    assert compare_exact_reals(lower, sp.sqrt(2)) < 0


def test_exact_truth_handles_algebraic_relation():
    assert exact_truth(sp.sqrt(2) > sp.Rational(7, 5)) is True
    assert exact_truth(sp.sqrt(2) < sp.Rational(7, 5)) is False


def test_root_isolation_declines_unsupported_domain():
    x = sp.Symbol("x")
    with pytest.raises(NotImplementedError, match="exact real-root isolation"):
        integral_real_roots(x**2 - sp.sqrt(2), x)


def test_string_resolution_rejects_ambiguous_same_name_symbols():
    x_plain = sp.Symbol("x")
    x_real = sp.Symbol("x", real=True)
    formula = sp.And(x_plain > 0, x_real > 0)
    with pytest.raises(ValueError, match="ambiguous"):
        semialg.simplify_boole(formula, ["x"])


def test_text_parser_reuses_explicit_variable_symbol_identity():
    from semialg.formula import parse_quant_form_text

    x = sp.Symbol("x", real=True)
    parsed = parse_quant_form_text("x > 0", variable_order=(x,))
    assert parsed.vars == (x,)
    assert parsed.matrix_expr.free_symbols == {x}


def test_text_parser_rejects_same_name_variable_ambiguity():
    from semialg.formula import parse_quant_form_text

    x_plain = sp.Symbol("x")
    x_real = sp.Symbol("x", real=True)
    with pytest.raises(ValueError, match="incompatible Symbol objects"):
        parse_quant_form_text("x > 0", variable_order=(x_plain, x_real))
