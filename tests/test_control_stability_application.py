import pytest
import sympy as sp

from semialg import equivalent
from semialg.applications import polynomial_stability_analysis, polynomial_stability_region


def test_quadratic_hurwitz_stability_region_is_exact():
    s, a, b = sp.symbols("s a b", real=True)
    result = polynomial_stability_analysis(s**2 + a * s + b, s, [a, b])
    assert result.certified
    assert result.degree == 2
    assert result.determinants == (a, a * b)
    assert equivalent(result.condition, sp.And(a > 0, b > 0), [a, b])


def test_stability_region_is_invariant_under_negative_scaling():
    s, a, b = sp.symbols("s a b", real=True)
    condition = polynomial_stability_region(-(s**2 + a * s + b), s, [a, b])
    assert equivalent(condition, sp.And(a > 0, b > 0), [a, b])


def test_cubic_hurwitz_determinants_match_classical_conditions():
    s, a0, a1, a2 = sp.symbols("s a0 a1 a2", real=True)
    result = polynomial_stability_analysis(s**3 + a2 * s**2 + a1 * s + a0, s)
    expected = (a2, a1 * a2 - a0, a0 * (a1 * a2 - a0))
    assert all(
        sp.expand(got - want) == 0 for got, want in zip(result.determinants, expected, strict=True)
    )
    assert result.parameters == (a0, a1, a2)


def test_stability_analysis_rejects_nonpolynomial_and_constant_inputs():
    s = sp.Symbol("s", real=True)
    with pytest.raises(ValueError, match="polynomial"):
        polynomial_stability_analysis(sp.sin(s) + 1, s)
    with pytest.raises(ValueError, match="positive degree"):
        polynomial_stability_analysis(sp.Integer(3), s)
