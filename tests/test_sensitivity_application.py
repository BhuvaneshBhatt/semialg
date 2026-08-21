import pytest
import sympy as sp

from semialg.applications import analyze_polynomial_sensitivity


def test_sensitivity_classifies_nondecreasing_and_constant_coordinates():
    x, y = sp.symbols("x y", real=True)
    domain = sp.And(x >= 0, x <= 2, y >= -1, y <= 1)
    result = analyze_polynomial_sensitivity(x**2 + 3 * y, [x, y], domain=domain)
    x_info = result.directions[x]
    y_info = result.directions[y]
    assert x_info.derivative == 2 * x
    assert x_info.classification == "nondecreasing"
    assert x_info.range_result.infimum == 0
    assert x_info.range_result.supremum == 4
    assert y_info.classification == "strictly_increasing"
    assert y_info.range_result.infimum == 3
    assert y_info.range_result.supremum == 3


def test_sensitivity_detects_mixed_derivative_sign():
    x = sp.Symbol("x", real=True)
    result = analyze_polynomial_sensitivity(x**2, [x], domain=sp.And(x >= -1, x <= 1))
    info = result.directions[x]
    assert info.classification == "mixed"
    assert not info.nonnegative
    assert not info.nonpositive


def test_sensitivity_detects_constant_model_coordinate():
    x, y = sp.symbols("x y", real=True)
    result = analyze_polynomial_sensitivity(
        x**2, [x, y], domain=sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    )
    assert result.directions[y].constant
    assert result.directions[y].classification == "constant"


def test_sensitivity_rejects_nonpolynomial_and_undeclared_parameters():
    x, beta = sp.symbols("x beta", real=True)
    with pytest.raises(ValueError, match="polynomial"):
        analyze_polynomial_sensitivity(sp.sin(x), [x])
    with pytest.raises(ValueError, match="undeclared symbolic parameters"):
        analyze_polynomial_sensitivity(beta * x + x**2, [x])
