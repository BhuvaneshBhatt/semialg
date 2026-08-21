import pytest
import sympy as sp

from semialg.applications import verify_polynomial_invariant


def test_inductive_safety_invariant_with_initial_and_unsafe_conditions():
    x = sp.Symbol("x", real=True)
    result = verify_polynomial_invariant(
        x >= 0,
        {x: x + 1},
        [x],
        initial_condition=sp.Eq(x, 0),
        unsafe_condition=x < 0,
    )
    assert result.certified
    assert result.inductive
    assert result.initial_valid is True
    assert result.safe_valid is True
    assert result.valid
    assert result.post_invariant == (x + 1 >= 0)


def test_noninductive_invariant_returns_counterexample():
    x = sp.Symbol("x", real=True)
    result = verify_polynomial_invariant(x >= 0, {x: x - 1}, [x])
    assert not result.inductive
    assert not result.valid
    assert result.counterexamples["inductive"] is not None


def test_sequence_transition_and_domain_are_supported():
    x, y = sp.symbols("x y", real=True)
    result = verify_polynomial_invariant(
        sp.And(x >= 0, y >= 0),
        (x + 1, y + x),
        [x, y],
        domain=sp.And(x >= 0, y >= 0),
    )
    assert result.inductive


def test_invariant_verifier_rejects_missing_or_nonpolynomial_updates():
    x, y = sp.symbols("x y", real=True)
    with pytest.raises(ValueError, match="missing"):
        verify_polynomial_invariant(x >= 0, {x: x + 1}, [x, y])
    with pytest.raises(ValueError, match="polynomial"):
        verify_polynomial_invariant(x >= 0, {x: sp.sin(x)}, [x])


def test_invariant_verifier_rejects_undeclared_symbolic_parameters():
    x, gain = sp.symbols("x gain", real=True)
    with pytest.raises(ValueError, match="undeclared symbolic parameters"):
        verify_polynomial_invariant(x >= 0, {x: gain * x}, [x])
