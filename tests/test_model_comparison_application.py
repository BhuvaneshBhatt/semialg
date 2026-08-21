import pytest
import sympy as sp

from semialg.applications import compare_polynomial_models


def test_model_comparison_reports_exact_discrepancy_and_dominance():
    x = sp.Symbol("x", real=True)
    domain = sp.And(x >= 0, x <= 1)
    result = compare_polynomial_models(x**2, x, [x], domain=domain)

    assert result.certified
    assert result.first_le_second
    assert not result.first_ge_second
    assert not result.equivalent_on_domain
    assert sp.simplify(result.minimum_difference.value + sp.Rational(1, 4)) == 0
    assert result.maximum_difference.value == 0
    assert sp.simplify(result.maximum_absolute_error - sp.Rational(1, 4)) == 0
    assert result.counterexamples["first_ge_second"] is not None


def test_model_comparison_certifies_algebraically_equivalent_models():
    x = sp.Symbol("x")
    result = compare_polynomial_models((x + 1) ** 2, x**2 + 2 * x + 1, ["x"], domain=x >= -2)
    assert result.variables == (x,)
    assert result.equivalent_on_domain
    assert result.maximum_absolute_error == 0


def test_model_comparison_rejects_nonpolynomial_and_undeclared_parameters():
    x, a = sp.symbols("x a", real=True)
    with pytest.raises(ValueError, match="polynomial"):
        compare_polynomial_models(sp.sin(x), x, [x])
    with pytest.raises(ValueError, match="undeclared"):
        compare_polynomial_models(a * x, x, [x])
