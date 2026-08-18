from __future__ import annotations

import sympy as sp

from semialg import (
    SignProofResult,
    prove_negative,
    prove_nonnegative,
    prove_nonpositive,
    prove_positive,
)


def test_prove_nonnegative_uses_square_certificate():
    x, y = sp.symbols("x y", real=True)
    result = prove_nonnegative(x**2 + y**2, [x, y], return_result=True)
    assert isinstance(result, SignProofResult)
    assert result.proven
    assert result.method == "sum_of_squares_or_even_powers"


def test_prove_positive_returns_counterexample_for_square():
    x = sp.symbols("x", real=True)
    result = prove_positive(x**2, [x], return_result=True)
    assert not result.proven
    assert result.counterexample == {x: 0}


def test_prove_positive_constant_plus_squares():
    x = sp.symbols("x", real=True)
    result = prove_positive(x**2 + 1, [x], return_result=True)
    assert result.proven
    assert result.method == "positive_constant_plus_squares"


def test_prove_nonpositive_negative_square():
    x = sp.symbols("x", real=True)
    result = prove_nonpositive(-(x**2), [x], return_result=True)
    assert result.proven


def test_boolean_api_preserved_for_sign_provers():
    x = sp.symbols("x", real=True)
    assert prove_nonnegative(x**2, [x]) is True
    assert prove_positive(x**2, [x]) is False
    assert prove_nonpositive(-(x**2), [x]) is True
    assert prove_negative(-(x**2), [x]) is False
