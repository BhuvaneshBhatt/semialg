from __future__ import annotations

import pytest
import sympy as sp

import semialg
from semialg.exact_arithmetic import compare_exact_reals, exact_truth
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
    result = semialg.normalize_domain_sensitive_constraints(sp.sqrt(x) >= 0, ["x"])
    assert result.variables == (x,)
    assert len(result.variables) == 1


def test_reasoning_string_variable_does_not_duplicate_assumption_variant():
    x = sp.Symbol("x")
    assert semialg.region_subset(x > 1, x > 0, ["x"])


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


def test_symbolic_integral_root_isolation_never_falls_back_to_nroots(monkeypatch):
    x = sp.Symbol("x")
    monkeypatch.setattr(
        sp,
        "real_roots",
        lambda *args, **kwargs: (_ for _ in ()).throw(NotImplementedError("forced")),
    )
    monkeypatch.setattr(
        sp,
        "nroots",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("nroots must not be used")),
    )
    with pytest.raises(NotImplementedError, match="exact real-root isolation"):
        integral_real_roots(x**2 - 2, x)


def test_string_resolution_rejects_ambiguous_same_name_symbols():
    x_plain = sp.Symbol("x")
    x_real = sp.Symbol("x", real=True)
    formula = sp.And(x_plain > 0, x_real > 0)
    with pytest.raises(ValueError, match="ambiguous"):
        semialg.simplify_boole(formula, ["x"])
