from __future__ import annotations

import pytest
import sympy as sp

from semialg import equivalent, implies, is_satisfiable

pytestmark = pytest.mark.contract

x, y = sp.symbols("x y", real=True)


@pytest.mark.parametrize(
    ("premise", "conclusion"),
    [
        (x > 1, x**2 > 1),
        (sp.And(x >= 0, x <= 1), x**2 <= 1),
        (x**2 + y**2 < 1, x**2 + y**2 <= 1),
        (sp.Eq(x, y), sp.Eq(x**2, y**2)),
    ],
)
def test_implication_agrees_with_counterexample_unsatisfiability(premise, conclusion) -> None:
    variables = tuple(
        sorted((premise.free_symbols | conclusion.free_symbols), key=lambda s: s.name)
    )
    expected = not is_satisfiable(sp.And(premise, sp.Not(conclusion)), variables)
    assert implies(premise, conclusion, variables) is expected


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (x**2 <= 1, sp.And(x >= -1, x <= 1)),
        (sp.Or(x < 0, x >= 0), sp.true),
        (sp.Eq((x - y) ** 2, 0), sp.Eq(x, y)),
    ],
)
def test_equivalence_agrees_with_two_implications(left, right) -> None:
    variables = tuple(sorted((left.free_symbols | right.free_symbols), key=lambda s: s.name))
    expected = implies(left, right, variables) and implies(right, left, variables)
    assert equivalent(left, right, variables) is expected


def test_decision_predicates_are_invariant_under_variable_renaming() -> None:
    u, v = sp.symbols("u v", real=True)
    formula = sp.And(x**2 + y**2 <= 4, x > y, y >= 0)
    renamed = formula.xreplace({x: u, y: v})
    assert is_satisfiable(formula, (x, y)) is is_satisfiable(renamed, (u, v))
    assert equivalent(formula, renamed.xreplace({u: x, v: y}), (x, y)) is True


def test_positive_polynomial_scaling_preserves_feasible_set() -> None:
    original = x**2 + y**2 <= 1
    scaled = 7 * (x**2 + y**2 - 1) <= 0
    assert equivalent(original, scaled, (x, y)) is True
