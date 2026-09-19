"""Neutral/sufficient/incompatible assumption triplets for public symbolic decisions."""

from __future__ import annotations

import pytest
import sympy as sp

import semialg


@pytest.mark.parametrize(
    ("assumptions", "expected"),
    (
        (sp.true, sp.Abs(sp.Symbol("x", real=True))),
        (sp.Symbol("x", real=True) >= 0, sp.Symbol("x", real=True)),
        (sp.Symbol("x", real=True) <= 0, -sp.Symbol("x", real=True)),
    ),
)
def test_simplification_assumption_matrix(assumptions, expected):
    x = sp.Symbol("x", real=True)
    assumptions = assumptions.xreplace({next(iter(assumptions.free_symbols), x): x})
    expected = expected.xreplace({next(iter(expected.free_symbols), x): x})
    assert semialg.simplify_under_assumptions(sp.sqrt(x**2), assumptions, (x,)) == expected


def test_real_valued_assumption_matrix_preserves_empty_domain_semantics():
    x = sp.Symbol("x", real=True)
    expression = sp.sqrt(x - 1)

    assert semialg.is_real_valued(expression, (x,), assumptions=sp.true) is False
    assert semialg.is_real_valued(expression, (x,), assumptions=x >= 1) is True
    assert semialg.is_real_valued(expression, (x,), assumptions=x < 1) is False
    assert semialg.is_real_valued(expression, (x,), assumptions=sp.false) is True


def test_function_sign_assumption_matrix_collapses_incompatible_domain():
    x = sp.Symbol("x", real=True)

    assert semialg.function_sign(x, (x,), assumptions=sp.true) == "mixed"
    assert semialg.function_sign(x, (x,), assumptions=x > 0) == "positive"
    assert semialg.function_sign(x, (x,), assumptions=x < 0) == "negative"
    assert semialg.function_sign(x, (x,), assumptions=sp.false) == "empty_domain"


def test_assumptions_keep_original_symbol_identity():
    x = sp.Symbol("x", real=True)
    result = semialg.simplify_under_assumptions(sp.sqrt(x**2), x >= 0, (x,))
    assert result == x
    assert result.free_symbols == {x}
    assert all(symbol.name == "x" and symbol.is_real is True for symbol in result.free_symbols)
