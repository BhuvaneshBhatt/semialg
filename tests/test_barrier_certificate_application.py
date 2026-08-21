import pytest
import sympy as sp

from semialg.applications import verify_barrier_certificate


def test_polynomial_barrier_certificate_is_certified():
    x = sp.Symbol("x", real=True)
    result = verify_barrier_certificate(
        x**2 - 1,
        {x: -x},
        [x],
        initial_condition=sp.Eq(x, 0),
        unsafe_condition=x**2 >= 4,
    )
    assert result.certified
    assert result.initial_valid
    assert result.unsafe_separated
    assert result.boundary_valid
    assert result.valid
    assert result.lie_derivative == -2 * x**2


def test_invalid_barrier_exposes_failed_separation():
    x = sp.Symbol("x", real=True)
    result = verify_barrier_certificate(
        x,
        {x: 0},
        [x],
        initial_condition=sp.Eq(x, 0),
        unsafe_condition=sp.Eq(x, 0),
    )
    assert not result.unsafe_separated
    assert not result.valid
    assert result.counterexamples["unsafe"] is not None


def test_strict_boundary_derivative_can_fail_when_nonstrict_passes():
    x = sp.Symbol("x", real=True)
    loose = verify_barrier_certificate(
        x**2 - 1,
        {x: 0},
        [x],
        initial_condition=sp.Eq(x, 0),
        unsafe_condition=x**2 >= 4,
    )
    strict = verify_barrier_certificate(
        x**2 - 1,
        {x: 0},
        [x],
        initial_condition=sp.Eq(x, 0),
        unsafe_condition=x**2 >= 4,
        derivative_strict=True,
    )
    assert loose.boundary_valid
    assert not strict.boundary_valid


def test_barrier_verifier_rejects_nonpolynomial_dynamics():
    x = sp.Symbol("x", real=True)
    with pytest.raises(ValueError, match="polynomial"):
        verify_barrier_certificate(
            x,
            {x: sp.sin(x)},
            [x],
            initial_condition=x <= 0,
            unsafe_condition=x >= 1,
        )
