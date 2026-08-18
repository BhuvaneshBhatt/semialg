from __future__ import annotations

import pytest
import sympy as sp

from semialg import (
    equivalent,
    implies,
    is_satisfiable,
    is_tautology,
    sample_point,
    sample_points,
    sign_at,
    sign_vector,
    solve_semialgebraic,
)
from semialg.decision import (
    EquivalenceResult,
    ImplicationResult,
    SatisfiabilityResult,
    TautologyResult,
)


def _assert_satisfies(formula: sp.Expr, point: dict[sp.Symbol, sp.Expr]) -> None:
    assert point is not None
    value = sp.simplify(formula.subs(point))
    assert value in (True, sp.true) or bool(value) is True


def test_decision_result_objects_match_boolean_api() -> None:
    x = sp.symbols("x", real=True)
    formula = (x >= 0) & (x <= 1)

    sat = is_satisfiable(formula, [x], return_result=True)
    assert isinstance(sat, SatisfiabilityResult)
    assert bool(sat) is is_satisfiable(formula, [x])
    assert sat.status == "sat"
    assert sat.witness is not None
    _assert_satisfies(formula, dict(sat.witness))

    taut = is_tautology(x**2 >= 0, [x], return_result=True)
    assert isinstance(taut, TautologyResult)
    assert bool(taut) is is_tautology(x**2 >= 0, [x])
    assert taut.counterexample is None

    implication = implies((x >= 0) & (x <= 1), x**2 <= 1, [x], return_result=True)
    assert isinstance(implication, ImplicationResult)
    assert bool(implication) is implies((x >= 0) & (x <= 1), x**2 <= 1, [x])
    assert implication.counterexample is None

    eq = equivalent(x**2 <= 1, (x >= -1) & (x <= 1), [x], return_result=True)
    assert isinstance(eq, EquivalenceResult)
    assert bool(eq) is equivalent(x**2 <= 1, (x >= -1) & (x <= 1), [x])
    assert eq.counterexample is None


def test_false_tautology_implication_and_equivalence_return_counterexamples() -> None:
    x = sp.symbols("x", real=True)

    taut = is_tautology(x**2 > 0, [x], return_result=True)
    assert not taut.valid
    assert taut.counterexample is not None
    _assert_satisfies(sp.Not(x**2 > 0), dict(taut.counterexample))

    implication = implies(x >= 0, x > 0, [x], return_result=True)
    assert not implication.valid
    assert implication.counterexample is not None
    _assert_satisfies((x >= 0) & sp.Not(x > 0), dict(implication.counterexample))

    eq = equivalent(x**2 < 1, (x >= -1) & (x <= 1), [x], return_result=True)
    assert not eq.equivalent
    assert eq.failed_direction in {"rhs_implies_lhs", "both"}
    assert eq.counterexample is not None


def test_solve_count_zero_preserves_satisfiability_without_sampling() -> None:
    x, y = sp.symbols("x y", real=True)
    triangle = (x >= 0) & (x <= 1) & (y >= x) & (y <= 1)

    result = solve_semialgebraic(triangle, [x, y], count=0)

    assert result.satisfiable
    assert result.samples == ()
    assert result.formula is not sp.false


def test_sampling_strategy_contracts_and_validation() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = (x >= 0) & (y >= 0) & (x + y <= 1)

    rational = sample_points(formula, [x, y], count=3, strategy="rational")
    assert len(rational) == 3
    assert all(point[x].is_Rational and point[y].is_Rational for point in rational)
    for point in rational:
        _assert_satisfies(formula, point)

    grid = sample_points(
        formula,
        [x, y],
        count=2,
        strategy="grid",
        bounds=[(0, 1), (0, 1)],
        grid_resolution=4,
    )
    assert len(grid) == 2
    for point in grid:
        _assert_satisfies(formula, point)

    random_left = sample_points(
        formula,
        [x, y],
        count=2,
        strategy="random",
        bounds=[(0, 1), (0, 1)],
        seed=7,
        exact=False,
    )
    random_right = sample_points(
        formula,
        [x, y],
        count=2,
        strategy="random",
        bounds=[(0, 1), (0, 1)],
        seed=7,
        exact=False,
    )
    assert random_left == random_right
    assert all(isinstance(point[x], sp.Float) for point in random_left)


def test_sample_point_uses_exact_algebraic_witness_when_available() -> None:
    x = sp.symbols("x", real=True)
    point = sample_point(sp.Eq(x**2, 2), [x], strategy="representative")

    assert point is not None
    assert sp.simplify(point[x] ** 2 - 2) == 0


def test_sign_api_exact_rational_algebraic_and_mapping_output() -> None:
    x = sp.symbols("x", real=True)

    assert sign_at(x**2 - 2, {x: sp.sqrt(2)}) == 0
    assert sign_at(x - 1, {x: sp.sqrt(2)}) == 1
    assert sign_at(x - 2, {x: sp.sqrt(2)}) == -1

    signs = sign_vector([x, x - 1, x**2 - 1], {x: 0}, as_dict=True)
    assert signs[x] == 0
    assert signs[x - 1] == -1
    assert signs[x**2 - 1] == -1


def test_unsupported_sampling_strategy_fails_clearly() -> None:
    x = sp.symbols("x", real=True)
    with pytest.raises(ValueError, match="unsupported sample strategy"):
        sample_points(x >= 0, [x], strategy="not-a-strategy")
