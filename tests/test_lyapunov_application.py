import pytest
import sympy as sp

from semialg.applications import verify_lyapunov_function


def test_quadratic_lyapunov_function_is_certified():
    x = sp.Symbol("x", real=True)
    result = verify_lyapunov_function(x**2, {x: -x}, [x])
    assert result.certified
    assert result.equilibrium_valid
    assert result.equilibrium_in_domain
    assert result.positive_definite
    assert result.derivative_valid
    assert result.valid
    assert result.lie_derivative == -2 * x**2


def test_bad_lyapunov_candidate_returns_counterexample():
    x = sp.Symbol("x", real=True)
    result = verify_lyapunov_function(x**2, {x: x}, [x])
    assert not result.derivative_valid
    assert not result.valid
    assert result.counterexamples["derivative"] is not None


def test_nonstrict_lyapunov_derivative_can_be_verified():
    x, y = sp.symbols("x y", real=True)
    result = verify_lyapunov_function(
        x**2 + y**2,
        {x: -x, y: 0},
        [x, y],
        derivative_strict=False,
    )
    assert result.positive_definite
    assert result.derivative_valid
    assert result.valid


def test_lyapunov_verifier_rejects_nonpolynomial_dynamics_and_parameters():
    x, gain = sp.symbols("x gain", real=True)
    with pytest.raises(ValueError, match="polynomial"):
        verify_lyapunov_function(x**2, {x: sp.sin(x)}, [x])
    with pytest.raises(ValueError, match="undeclared symbolic parameters"):
        verify_lyapunov_function(x**2, {x: -gain * x}, [x])


def test_lyapunov_equilibrium_must_lie_in_domain():
    x = sp.Symbol("x", real=True)
    result = verify_lyapunov_function(x**2, {x: -x}, [x], domain=x >= 1)
    assert not result.equilibrium_in_domain
    assert not result.equilibrium_valid
    assert not result.valid
