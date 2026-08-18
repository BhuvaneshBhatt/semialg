import pytest
import sympy as sp

from semialg import integrate_over_region

pytestmark = pytest.mark.slow


def test_integrate_over_region_symbolic_is_default():
    x = sp.symbols("x", real=True)
    assert sp.simplify(integrate_over_region(x**2, x**2 <= 1, [x]) - sp.Rational(2, 3)) == 0
    assert (
        sp.simplify(
            integrate_over_region(x**2, x**2 <= 1, [x], method="symbolic") - sp.Rational(2, 3)
        )
        == 0
    )


def test_integrate_over_region_numeric_method_returns_approximate_result():
    x = sp.symbols("x", real=True)
    result = integrate_over_region(
        sp.cos(x**3 + x), (x >= 0) & (x <= 1), [x], method="numeric", return_result=True
    )
    expected = sp.Integral(sp.cos(x**3 + x), (x, 0, 1)).evalf(50)
    assert result.exact is False
    assert result.diagnostics["evaluation_method"] == "numeric"
    assert abs(result.value - expected) < sp.Float("1e-40")


def test_integrate_over_region_auto_falls_back_to_numeric_when_symbolic_fails():
    x = sp.symbols("x", real=True)
    with pytest.raises(NotImplementedError):
        integrate_over_region(sp.cos(x**3 + x), (x >= 0) & (x <= 1), [x], method="symbolic")
    result = integrate_over_region(
        sp.cos(x**3 + x), (x >= 0) & (x <= 1), [x], method="auto", return_result=True
    )
    assert result.exact is False
    assert result.diagnostics["evaluation_method"] == "numeric"
    assert 0 < result.value < 1


def test_integrate_over_region_rejects_unknown_method():
    x = sp.symbols("x", real=True)
    with pytest.raises(ValueError):
        integrate_over_region(x, (x >= 0) & (x <= 1), [x], method="hybrid")
