from __future__ import annotations

import sympy as sp

from semialg.preprocess import semialgebraicize


def test_odd_denominator_power_does_not_require_nonnegative_base():
    x = sp.symbols("x", real=True)
    result = semialgebraicize(x ** sp.Rational(1, 3) <= 2, variables=(x,))
    assert result.aux_vars
    aux = result.aux_vars[0]
    assert sp.Eq(aux**3, x) in result.assumptions
    assert not any(cond == (x >= 0) for cond in result.assumptions)


def test_even_numerator_odd_denominator_power_records_nonnegative_value():
    x = sp.symbols("x", real=True)
    result = semialgebraicize(x ** sp.Rational(2, 3) <= 2, variables=(x,))
    aux = result.aux_vars[0]
    assert (aux >= 0) in result.assumptions
    assert sp.Eq(aux**3, x**2) in result.assumptions


def test_abs_fractional_power_creates_abs_and_power_auxiliaries():
    x, y = sp.symbols("x y", real=True)
    result = semialgebraicize(
        sp.Abs(x) ** sp.Rational(3, 2) + sp.Abs(y) ** sp.Rational(3, 2) <= 1, variables=(x, y)
    )
    kinds = [aux.kind for aux in result.auxiliary_defs]
    assert kinds.count("abs") == 2
    assert kinds.count("rational_power") == 2
