from __future__ import annotations

import sympy as sp

from semialg.solve.integer.congruence import norm_mod_form, solve_quant_free_mod_sys


def test_normalized_mod_residual_uses_prime_linear_fast_path():
    x, y = sp.symbols("x y", integer=True)
    formula = sp.And(sp.Eq(x + y, 1), sp.Eq(x - y, 3))

    normalized = norm_mod_form(formula, 5)
    result = solve_quant_free_mod_sys(normalized, (x, y), 5)

    assert result.complete
    assert result.method == "linear_row_reduction_mod_prime"
    assert result.points == [(2, 4)]


def test_normalized_mod_residual_uses_composite_linear_fast_path():
    x, y = sp.symbols("x y", integer=True)
    formula = sp.And(sp.Eq(x + y, 1), sp.Eq(x - y, 3))

    normalized = norm_mod_form(formula, 8)
    result = solve_quant_free_mod_sys(normalized, (x, y), 8)

    assert result.complete
    assert result.method == "smith_normal_form_mod_composite"
    assert set(result.points) == {(2, 7), (6, 3)}
