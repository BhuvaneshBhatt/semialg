from __future__ import annotations

import itertools

import pytest
import sympy as sp

from semialg import solvability_conditions
from semialg.parameters import root_count_conditions


def _selected_count(conditions, assignment):
    selected = [
        count
        for count, condition in conditions.items()
        if sp.simplify(condition.subs(assignment)) in (sp.true, True)
    ]
    assert len(selected) == 1
    return selected[0]


def _distinct_real_root_count(expr, variable):
    return len(tuple(dict.fromkeys(sp.real_roots(sp.Poly(expr, variable)))))


@pytest.mark.parametrize(
    "coefficients",
    [
        (1, 0, -5, 0, 4),  # four real roots
        (1, 0, 1, 0, 1),  # no real roots
        (1, 0, 0, 0, -1),  # two real roots
        (2, 1, -2, -1, 1),
        (-1, 2, 1, -2, 3),
        (0, 1, -3, 2, 0),  # cubic degree drop
        (0, 0, 1, -2, 1),  # quadratic degree drop
        (0, 0, 0, 2, -4),  # linear degree drop
        (0, 0, 0, 0, 5),  # nonzero constant
    ],
)
def test_quartic_certified_strata_match_exact_specialization(coefficients):
    x, a, b, c, d, e = sp.symbols("x a b c d e", real=True)
    family = a * x**4 + b * x**3 + c * x**2 + d * x + e
    conditions = root_count_conditions(family, x, (a, b, c, d, e))
    assignment = dict(zip((a, b, c, d, e), coefficients, strict=True))
    specialized = sp.expand(family.subs(assignment))
    selected = _selected_count(conditions, assignment)

    if selected == -1:
        assert sp.discriminant(specialized, x) == 0
    elif specialized == 0:
        assert selected is sp.oo
    else:
        assert selected == _distinct_real_root_count(specialized, x)


def test_quartic_positive_discriminant_root_counts():
    x, a = sp.symbols("x a", real=True)
    family = x**4 + a * x**2 + 1
    conditions = root_count_conditions(family, x, (a,))

    assert _selected_count(conditions, {a: -3}) == 4
    assert _selected_count(conditions, {a: 0}) == 0
    assert _selected_count(conditions, {a: 3}) == 0


def test_quartic_multiple_root_locus_remains_explicitly_unknown():
    x, a = sp.symbols("x a", real=True)
    family = x**4 + a * x**2 + 1
    conditions = root_count_conditions(family, x, (a,))

    assert _selected_count(conditions, {a: -2}) == -1
    assert _selected_count(conditions, {a: 2}) == -1


def test_partial_quartic_stratified_result_does_not_certify_unknown_locus():
    x, a = sp.symbols("x a", real=True)
    family = x**4 + a * x**2 + 1

    result = root_count_conditions(family, x, (a,), return_stratified=True)

    assert result.complete is False
    assert sp.simplify(result.diagnostics["unknown_condition"].subs(a, -2)) is sp.true
    assert all(branch.value != -1 for branch in result.branches)


def test_quartic_unknown_root_count_falls_back_for_existence_query():
    x, a = sp.symbols("x a", real=True)
    family = x**4 + a * x**2 + 1

    condition = solvability_conditions(sp.Eq(family, 0), (x,), (a,))

    assert sp.simplify(condition.subs(a, -3)) is sp.true
    assert sp.simplify(condition.subs(a, -2)) is sp.true
    assert sp.simplify(condition.subs(a, 0)) is sp.false
    assert sp.simplify(condition.subs(a, 3)) is sp.false


def test_small_integer_quartics_never_misclassify_a_certified_stratum():
    x = sp.symbols("x", real=True)
    a, b, c, d, e = sp.symbols("a b c d e", real=True)
    family = a * x**4 + b * x**3 + c * x**2 + d * x + e
    conditions = root_count_conditions(family, x, (a, b, c, d, e))

    # Exhaust a compact but structurally varied coefficient grid. Unknown
    # discriminant-zero strata are permitted; certified counts are not.
    values = (-1, 0, 1)
    for coefficients in itertools.product(values, repeat=5):
        if coefficients == (0, 0, 0, 0, 0):
            continue
        assignment = dict(zip((a, b, c, d, e), coefficients, strict=True))
        selected = _selected_count(conditions, assignment)
        if selected == -1:
            continue
        specialized = sp.expand(family.subs(assignment))
        if specialized == 0:
            assert selected is sp.oo
        else:
            assert selected == _distinct_real_root_count(specialized, x)
