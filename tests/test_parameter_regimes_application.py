import pytest
import sympy as sp

from semialg.applications import analyze_parameter_regimes, analyze_root_count_regimes


def test_parameter_regimes_partition_solvability_exactly():
    x, a = sp.symbols("x a", real=True)
    result = analyze_parameter_regimes(sp.Eq(x**2 + a, 0), [x], [a])

    assert result.certified
    assert result.quantity == "solvability"
    assert result.select({a: -1}) is True
    assert result.select({a: 1}) is False
    assert result.regime_count == 2


def test_root_count_regimes_capture_discriminant_changes():
    x, a = sp.symbols("x a", real=True)
    result = analyze_root_count_regimes(x**2 + a, x, [a])

    assert result.certified
    assert result.quantity == "real_root_count"
    assert result.select({a: -1}) == 2
    assert result.select({a: 0}) == 1
    assert result.select({a: 1}) == 0
    assert result.regime_count == 3


def test_parameter_regimes_preserve_string_symbol_identity_and_validate_symbols():
    x = sp.Symbol("x")
    a = sp.Symbol("a")
    result = analyze_parameter_regimes(sp.Eq(x, a), ["x"], ["a"])
    assert result.variables == (x,)
    assert result.parameters == (a,)

    with pytest.raises(ValueError, match="disjoint"):
        analyze_parameter_regimes(x >= 0, [x], [x])

    b = sp.Symbol("b")
    with pytest.raises(ValueError, match="declared"):
        analyze_parameter_regimes(sp.Eq(x, a + b), [x], [a])


def test_root_count_regimes_reject_nonpolynomial_input():
    x, a = sp.symbols("x a", real=True)
    with pytest.raises(ValueError, match="polynomial"):
        analyze_root_count_regimes(sp.sin(x) + a, x, [a])
