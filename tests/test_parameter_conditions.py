from __future__ import annotations

import sympy as sp

from semialg import root_count_conditions, solvability_conditions


def test_solvability_conditions_for_quadratic_real_root_parameters() -> None:
    x, a, b = sp.symbols("x a b", real=True)
    condition = solvability_conditions(sp.Eq(x**2 + a * x + b, 0), [x], [a, b])
    assert sp.simplify(condition) == (a**2 - 4 * b >= 0)


def test_solvability_conditions_for_strict_negative_quadratic_parameter() -> None:
    x, a = sp.symbols("x a", real=True)
    condition = solvability_conditions(x**2 + a < 0, [x], [a])
    assert sp.simplify(condition) == (a < 0)


def test_solvability_conditions_can_return_result_object() -> None:
    x, a = sp.symbols("x a", real=True)
    result = solvability_conditions(sp.Eq(x**2, a), [x], [a], return_result=True)
    assert result.formula == (a >= 0)
    assert result.variables == (x,)
    assert result.parameters == (a,)
    assert bool(result) is True


def test_root_count_conditions_for_monic_quadratic_family() -> None:
    x, a, b = sp.symbols("x a b", real=True)
    conditions = root_count_conditions(x**2 + a * x + b, x, [a, b])
    assert conditions[2] == (a**2 - 4 * b > 0)
    assert conditions[1] == sp.Eq(a**2 - 4 * b, 0)
    assert conditions[0] == (a**2 - 4 * b < 0)


def test_root_count_conditions_result_preserves_classification() -> None:
    x, a = sp.symbols("x a", real=True)
    result = root_count_conditions(x**2 + a, x, [a], return_result=True)
    assert result.parameters == (a,)
    assert result.condition_for_count(2) == (a < 0)
    assert result.condition_for_count(1) == sp.Eq(a, 0)
    assert result.condition_for_count(0) == (a > 0)
    assert result.classification.cells


def test_solvability_string_names_reuse_existing_unassumed_symbols() -> None:
    x, a = sp.symbols("x a")
    result = solvability_conditions(sp.Eq(x**2, a), ["x"], ["a"], return_stratified=True)
    assert result.parameters == (a,)
    assert result.select({"a": 4}) is True
    assert result.select({"a": -1}) is False
    assert all(x not in branch.condition.free_symbols for branch in result.branches)
