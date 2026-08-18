from __future__ import annotations

import sympy as sp

from semialg import is_satisfiable, sample_point, solve_semialgebraic


def test_solve_semialgebraic_auto_uses_rur_for_zero_dimensional_equalities():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x**2 + y**2, 1), sp.Eq(x, y))

    result = solve_semialgebraic(formula, [x, y], count=10)

    assert result.method == "rational_univariate"
    assert result.satisfiable
    assert len(result.samples) == 2
    assert all(sp.simplify(point[x] - point[y]) == 0 for point in result.samples)
    assert result.diagnostics["used_rur"] is True


def test_solve_semialgebraic_rur_filters_inequality_branch_exactly():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x**2 + y**2, 1), sp.Eq(x, y), x > 0)

    result = solve_semialgebraic(formula, [x, y], method="rur", count=10)

    assert result.method == "rational_univariate"
    assert result.satisfiable
    assert len(result.samples) == 1
    assert sp.simplify(result.samples[0][x] - sp.sqrt(2) / 2) == 0


def test_is_satisfiable_returns_rur_witness_for_finite_algebraic_system():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x**2 + y**2, 1), sp.Eq(x, y), x > 0)

    result = is_satisfiable(formula, [x, y], return_result=True)

    assert result.satisfiable
    assert result.method == "rational_univariate"
    assert result.witness is not None
    assert sp.simplify(result.witness[x] - sp.sqrt(2) / 2) == 0


def test_sample_point_uses_rur_for_irrational_finite_witness():
    x = sp.symbols("x", real=True)
    point = sample_point(sp.Eq(x**2, 2), [x])

    assert point is not None
    assert sp.simplify(point[x] ** 2 - 2) == 0
