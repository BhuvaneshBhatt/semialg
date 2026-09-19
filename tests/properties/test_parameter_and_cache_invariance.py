from __future__ import annotations

import itertools

import pytest
import sympy as sp
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg import is_satisfiable, solve_semialgebraic
from semialg.algebraic.cache import algebraic_cache_stats, clear_algebraic_caches
from semialg.algebraic.comparison import compare_samples
from semialg.algebraic.intervals import RationalInterval
from semialg.algebraic.roots import isolate_real_roots, refine_isol_intvl
from semialg.algebraic.samples import AlgebraicRoot, RationalSample

X = sp.Symbol("x", real=True)
A, B, C = sp.symbols("a b c", real=True)


@pytest.fixture(scope="module")
def linear_condition():
    return solve_semialgebraic(
        sp.Eq(A * X + B, 0),
        (X,),
        parameters=(A, B),
        count=0,
        output="conditions",
    )


@pytest.fixture(scope="module")
def quadratic_condition():
    return solve_semialgebraic(
        sp.Eq(A * X**2 + B * X + C, 0),
        (X,),
        parameters=(A, B, C),
        count=0,
        output="conditions",
    )


def _truth(value: sp.Expr) -> bool:
    simplified = sp.simplify(value)
    assert simplified in (sp.true, sp.false)
    return simplified is sp.true


@settings(max_examples=50, deadline=None)
@given(a=st.integers(-4, 4), b=st.integers(-4, 4))
def test_fast_linear_parameter_condition_matches_specialized_decision(
    a: int, b: int, linear_condition
) -> None:
    predicted = _truth(linear_condition.subs({A: a, B: b}))
    observed = is_satisfiable(sp.Eq(a * X + b, 0), (X,))
    assert predicted is observed


@settings(max_examples=60, deadline=None)
@given(
    a=st.integers(-3, 3),
    b=st.integers(-3, 3),
    c=st.integers(-3, 3),
)
def test_fast_quadratic_parameter_condition_matches_specialized_decision(
    a: int, b: int, c: int, quadratic_condition
) -> None:
    predicted = _truth(quadratic_condition.subs({A: a, B: b, C: c}))
    observed = is_satisfiable(sp.Eq(a * X**2 + b * X + c, 0), (X,))
    assert predicted is observed


@pytest.mark.parametrize("use_strings", [False, True])
@pytest.mark.parametrize("order", [("x", "y"), ("y", "x")])
def test_solve_symbol_resolution_is_invariant_across_names_and_orders(
    use_strings: bool, order: tuple[str, str]
) -> None:
    x = sp.Symbol("x")
    y = sp.Symbol("y", real=True)
    a = sp.Symbol("a", nonnegative=True)
    by_name = {symbol.name: symbol for symbol in (x, y)}
    variables = order if use_strings else tuple(by_name[name] for name in order)
    parameters = ("a",) if use_strings else (a,)
    result = solve_semialgebraic(
        sp.And(sp.Eq(x + y, a), x >= 0),
        variables,
        parameters=parameters,
        variable_order=variables,
        count=0,
    )
    assert result.variables == tuple(by_name[name] for name in order)
    assert result.parameters == (a,)


def test_all_solve_order_hints_reject_ambiguous_string_symbols() -> None:
    x_plain = sp.Symbol("x")
    x_real = sp.Symbol("x", real=True)
    formula = sp.And(x_plain >= 0, x_real <= 1)
    for keyword in ("variable_order", "projection_order"):
        with pytest.raises(ValueError, match="ambiguous"):
            solve_semialgebraic(formula, (x_plain, x_real), **{keyword: ("x",)}, count=0)


def test_comparison_cache_is_independent_of_history_and_operand_orientation() -> None:
    x = sp.Symbol("x", real=True)
    polynomial = sp.Poly(x**2 - 2, x)
    positive = AlgebraicRoot(polynomial, RationalInterval(1, 2), root_index=1)
    negative = AlgebraicRoot(polynomial, RationalInterval(-2, -1), root_index=0)
    zero = RationalSample(0)
    expected = {(id(negative), id(zero)): -1, (id(zero), id(positive)): -1}

    comparisons = ((negative, zero), (zero, positive))
    for history in itertools.permutations(comparisons):
        clear_algebraic_caches()
        for left, right in history:
            assert compare_samples(left, right) == expected[(id(left), id(right))]
            assert compare_samples(right, left) == -expected[(id(left), id(right))]


def test_refined_intervals_retain_semantic_comparison_cache_reuse() -> None:
    x = sp.Symbol("x", real=True)
    roots = isolate_real_roots(sp.Poly(x**2 - 2, x))
    left, right = roots
    refined_left = refine_isol_intvl(left, steps=2)
    refined_right = refine_isol_intvl(right, steps=2)

    clear_algebraic_caches()
    assert compare_samples(left, right) == -1
    before = algebraic_cache_stats()
    assert compare_samples(refined_left, refined_right) == -1
    after = algebraic_cache_stats()
    # Disjoint certificates take the structural interval path. Refinement may
    # miss that cheap cache, but must never require algebraic comparison work.
    assert after.comparison_refinements == before.comparison_refinements == 0
