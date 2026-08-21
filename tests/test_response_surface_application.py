import pytest
import sympy as sp

from semialg import equivalent
from semialg.applications import analyze_response_surface


def test_polynomial_response_surface_analysis_on_box():
    x, y, t = sp.symbols("x y t", real=True)
    domain = sp.And(x >= -1, x <= 1, y >= -1, y <= 1)
    result = analyze_response_surface(x**2 + y**2, [x, y], domain=domain, thresholds=[1])
    assert result.certified
    assert result.minimum.value == 0
    assert result.maximum.value == 2
    assert result.gradient == (2 * x, 2 * y)
    assert result.stationary_condition == sp.And(sp.Eq(2 * x, 0), sp.Eq(2 * y, 0))
    assert equivalent(result.range_result.formula, sp.And(t >= 0, t <= 2), [t])
    assert result.threshold_regions[sp.Integer(1)] == sp.And(domain, x**2 + y**2 >= 1)


def test_response_surface_preserves_unassumed_symbol_identity():
    x = sp.Symbol("x")
    result = analyze_response_surface(x**2, ["x"], domain=sp.And(x >= -1, x <= 1))
    assert result.variables == (x,)
    assert result.minimum.value == 0
    assert result.maximum.value == 1


def test_response_surface_rejects_nonpolynomial_models_and_empty_variables():
    x = sp.Symbol("x", real=True)
    with pytest.raises(ValueError, match="polynomial"):
        analyze_response_surface(sp.sin(x), [x])
    with pytest.raises(ValueError, match="predictor"):
        analyze_response_surface(sp.Integer(1), [])


def test_response_surface_rejects_undeclared_symbolic_parameters():
    x, beta = sp.symbols("x beta", real=True)
    with pytest.raises(ValueError, match="undeclared symbolic parameters"):
        analyze_response_surface(beta * x + x**2, [x], domain=sp.And(x >= -1, x <= 1))
