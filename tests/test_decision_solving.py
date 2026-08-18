from __future__ import annotations

import pytest
import sympy as sp

from semialg import (
    equivalent,
    implies,
    is_satisfiable,
    is_tautology,
    solve_semialgebraic,
)


def test_is_satisfiable_detects_basic_real_sets():
    x = sp.Symbol("x", real=True)

    assert is_satisfiable(x**2 <= 1, [x]) is True
    assert is_satisfiable(sp.And(x > 0, x < 0), [x]) is False
    assert is_satisfiable(sp.Eq(x**2 + 1, 0), [x]) is False


def test_is_tautology_uses_real_semantics():
    x = sp.Symbol("x", real=True)

    assert is_tautology(x**2 >= 0, [x]) is True
    assert is_tautology(x**2 > 0, [x]) is False
    assert is_tautology(sp.Or(x <= 0, x > 0), [x]) is True


def test_implication_with_assumptions():
    x, y = sp.symbols("x y", real=True)

    assert implies(x > 1, x**2 > 1, [x]) is True
    assert implies(sp.And(x >= 0, y >= 0), x * y >= 0, [x, y]) is True
    assert implies(x >= 0, x > 0, [x]) is False


def test_equivalent_rewrites_common_interval_formula():
    x = sp.Symbol("x", real=True)

    assert equivalent(x**2 <= 1, sp.And(x >= -1, x <= 1), [x]) is True
    assert equivalent(x**2 < 1, sp.And(x >= -1, x <= 1), [x]) is False


def test_solve_semialgebraic_reduces_univariate_interval():
    x = sp.Symbol("x", real=True)

    result = solve_semialgebraic(x**2 <= 1, [x])

    assert result.satisfiable is True
    assert result.formula == sp.And(x >= -1, x <= 1)
    assert result.sample is not None
    assert result.sample[x] == 0


def test_solve_semialgebraic_reports_unsatisfiable_system():
    x = sp.Symbol("x", real=True)

    result = solve_semialgebraic([x > 0, x < 0], [x])

    assert result.satisfiable is False
    assert result.formula is sp.false
    assert result.sample is None
    assert result.samples == ()


@pytest.mark.slow
def test_solve_semialgebraic_accepts_constraint_lists_and_returns_samples():
    x, y = sp.symbols("x y", real=True)

    result = solve_semialgebraic([x >= 0, y >= 0, x + y <= 1], [x, y], count=2)

    assert result.satisfiable is True
    assert result.variables == (x, y)
    assert len(result.samples) >= 1
    for sample in result.samples:
        assert sample[x] >= 0
        assert sample[y] >= 0
        assert sample[x] + sample[y] <= 1
