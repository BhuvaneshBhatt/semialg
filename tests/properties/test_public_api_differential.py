"""Bounded generated problems checked against exact independent oracles."""

from __future__ import annotations

from collections import Counter

import pytest
import sympy as sp
from hypothesis import example, given, settings
from hypothesis import strategies as st

from semialg import cad, classify_real_roots, is_satisfiable, solve_semialgebraic
from semialg.algebraic.cache import clear_algebraic_caches
from semialg.cad_algorithms.performance_cache import clear_cad_caches
from semialg.formula import parse_formula
from semialg.qe.complete import qe_by_complete_cad

pytestmark = pytest.mark.slow

X, Y, U, V = sp.symbols("x y u v", real=True)
A, B, C = sp.symbols("a b c", real=True)
SMALL = st.integers(-3, 3)


@pytest.fixture(scope="module")
def quadratic_family():
    return classify_real_roots(A * X**2 + B * X + C, X, parameters=(A, B, C))


@settings(max_examples=24, deadline=None)
@given(a=SMALL, b=SMALL, c=SMALL)
@example(a=0, b=0, c=0)
@example(a=0, b=0, c=1)
@example(a=0, b=1, c=0)
@example(a=1, b=2, c=1)
def test_parameter_root_counts_match_specialization_and_discriminant(quadratic_family, a, b, c):
    discriminant = b * b - 4 * a * c
    expected = (
        (2 if discriminant > 0 else 1 if discriminant == 0 else 0)
        if a
        else (1 if b else 0 if c else sp.oo)
    )
    selected = [
        cell
        for cell in quadratic_family.cells
        if sp.simplify(cell.condition.subs({A: a, B: b, C: c})) is sp.true
    ]
    assert len(selected) == 1
    assert selected[0].root_count == expected
    specialized = classify_real_roots(a * X**2 + b * X + c, X)
    assert specialized.generic_root_count == expected


@settings(max_examples=24, deadline=None)
@given(roots=st.lists(SMALL, min_size=1, max_size=4), unit=st.sampled_from((-2, -1, 1, 2)))
@example(roots=[0, 0, 0, 0], unit=-1)
def test_root_classification_matches_generated_factor_multiplicities(roots, unit):
    polynomial = unit * sp.prod(X - root for root in roots)
    expected = tuple(sorted(Counter(roots).values()))
    original = classify_real_roots(polynomial, X)
    renamed = classify_real_roots(sp.expand(polynomial).xreplace({X: U}), U)
    for result in (original, renamed):
        assert result.generic_root_count == len(set(roots))
        assert result.generic_multiplicity_pattern == expected


@settings(max_examples=20, deadline=None)
@given(left=SMALL, width=st.integers(0, 4), strict=st.booleans(), multiplier=st.integers(1, 4))
@example(left=0, width=0, strict=True, multiplier=1)
@example(left=0, width=0, strict=False, multiplier=2)
def test_solver_and_cad_match_factored_quadratic_sets(left, width, strict, multiplier):
    polynomial = multiplier * (X - left) * (X - left - width)
    formula = polynomial < 0 if strict else polynomial <= 0
    expected = sp.Interval(left, left + width, strict, strict)
    solved = solve_semialgebraic(formula, (X,), count=0)
    decomposed = cad(formula, (X,), return_result=True)
    assert decomposed.status == "complete"
    assert solved.formula.as_set() & sp.S.Reals == expected
    assert decomposed.formula.as_set() & sp.S.Reals == expected
    assert is_satisfiable(formula, (X,)) is (expected != sp.S.EmptySet)


@settings(max_examples=16, deadline=None)
@given(a=SMALL, b=SMALL, lower=SMALL, upper=SMALL)
@example(a=0, b=0, lower=0, upper=0)
@example(a=2, b=-1, lower=1, upper=0)
def test_affine_feasibility_matches_elimination_under_renaming_and_order(a, b, lower, upper):
    formula = sp.And(sp.Eq(Y, a * X + b), X >= lower, X <= upper)
    expected = lower <= upper
    variants = ((formula, (X, Y)), (formula, (Y, X)), (formula.xreplace({X: U, Y: V}), (V, U)))
    for condition, variables in variants:
        decision = is_satisfiable(condition, variables, return_result=True)
        assert decision.satisfiable is expected
        if expected:
            assert decision.witness is not None
            assert sp.simplify(condition.subs(decision.witness)) is sp.true
    # Disable presolve and the variety shortcut for the CAD differential path.
    eliminated = qe_by_complete_cad(
        (X, Y),
        (("exists", X), ("exists", Y)),
        parse_formula(formula),
        free_variables=(),
        variable_order_strategy="preserve",
        use_presolve=False,
        allow_variety_cad=False,
        return_result=True,
    )
    assert eliminated.truth_value is expected


@settings(max_examples=16, deadline=None)
@given(center=SMALL, radius=st.integers(0, 3), strict=st.booleans())
@example(center=0, radius=0, strict=True)
@example(center=0, radius=0, strict=False)
def test_decisions_survive_cold_warm_and_unrelated_cache_history(center, radius, strict):
    polynomial = (X - center) ** 2 - radius**2
    formula = polynomial < 0 if strict else polynomial <= 0
    expected = radius > 0 or not strict
    clear_algebraic_caches()
    clear_cad_caches()
    try:
        cold = is_satisfiable(formula, (X,))
        warm = is_satisfiable(formula, (X,))
        assert is_satisfiable(sp.Eq(X**2, -1), (X,)) is False
        assert is_satisfiable(sp.Eq(U**2, 2), (U,)) is True
        polluted = is_satisfiable(formula, (X,))
        clear_algebraic_caches()
        clear_cad_caches()
        fresh = is_satisfiable(formula, (X,))
        assert cold is warm is polluted is fresh is expected
    finally:
        clear_algebraic_caches()
        clear_cad_caches()
