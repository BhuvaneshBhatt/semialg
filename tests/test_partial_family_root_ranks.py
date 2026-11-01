"""Partial root families must not become complete boundary-rank certificates."""

import pytest
import sympy as sp

from semialg.algebraic.cache import clear_algebraic_caches
from semialg.algebraic.comparison import (
    _rank_from_known_root_intervals,
    compare_samples,
)
from semialg.algebraic.intervals import RationalInterval
from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext, RationalSample
from semialg.cad_algorithms.decomposition import _stack_roots_over_point
from semialg.context import computation_context


@pytest.fixture(autouse=True)
def exact_context():
    with computation_context():
        yield


def _selected_roots():
    x, y = sp.symbols("x y", real=True)
    polynomial = (y - 1) * (y - 2) * (y - 3)
    left = AlgebraicRoot(
        sp.Poly(polynomial, y),
        RationalInterval(sp.Rational(3, 2), sp.Rational(5, 2)),
        1,
        1,
        None,
        FiberRootContext(sp.Poly(polynomial + x, x, y), (x, y), (RationalSample(sp.S.Zero),), y),
    )
    right = AlgebraicRoot(
        sp.Poly(10 * y - 21, y),
        RationalInterval(sp.S(2), sp.S(8)),
        0,
        1,
        None,
        FiberRootContext(sp.Poly(10 * y - 21 + x, x, y), (x, y), (RationalSample(sp.S.Zero),), y),
    )
    return left, right


@pytest.mark.parametrize("reverse", [False, True])
def test_selected_root_order(reverse):
    clear_algebraic_caches()
    left, right = _selected_roots()
    pair = (right, left) if reverse else (left, right)
    expected = 1 if reverse else -1  # 2 < 21/10, with valid overlapping intervals.
    assert compare_samples(*pair) == expected
    assert compare_samples(*pair) == expected
    # A comparison must not poison the family's exact counts for later work.
    assert left.fiber_context.roots_left_of(4) == 3


def test_partial_nodes_do_not_count_siblings():
    left, _ = _selected_roots()
    assert _rank_from_known_root_intervals((left,), sp.S(3)) is None
    assert _rank_from_known_root_intervals((left,), sp.S.Zero) is None
    assert _rank_from_known_root_intervals((), sp.S.Zero) is None


def test_adjacent_absolute_ranks_determine_boundary():
    left, _ = _selected_roots()
    upper = AlgebraicRoot(left.polynomial, RationalInterval(sp.S(3), sp.S(3)), 2)
    assert _rank_from_known_root_intervals((left, upper), sp.Rational(11, 4)) == 2
    assert _rank_from_known_root_intervals((left, upper), sp.S(4), total_roots=3) == 3


def test_family_order_survives_unrelated_merges():
    clear_algebraic_caches()
    x, y = sp.symbols("x y", real=True)
    left = 27 * x**4 + 64 * x**3 + 192 * x**2 * y + 192 * x * y**2 + 64 * y**3
    right = (
        4 * x**5
        + 13 * x**4
        + 104 * x**3 * y
        + 32 * x**3
        + 16 * x**2 * y**2
        + 224 * x**2 * y
        + 640 * x * y**2
        + 512 * x * y
        + 256 * y**3
        + 256 * y**2
    )
    prefix = (RationalSample(sp.S(-9)),)
    variables = (x, y)
    families = (sp.Poly(left, *variables), sp.Poly(right, *variables))
    first = _stack_roots_over_point(families[:1], variables, prefix, y)[0]
    second = _stack_roots_over_point(families[1:], variables, prefix, y)[1]
    separator = -sp.Rational(101, 20)
    assert sp.Poly(left.subs(x, -9), y).count_roots(-sp.oo, separator) == 0
    assert sp.Poly(right.subs(x, -9), y).count_roots(-sp.oo, separator) == 2
    _stack_roots_over_point(families, variables, prefix, y)
    assert compare_samples(first, second) == 1


def test_ranked_node_excludes_previous_endpoint():
    clear_algebraic_caches()
    x, y = sp.symbols("x y", real=True)
    polynomial = (y - 1) * (y - 2) * (y - 3)
    context = FiberRootContext(
        sp.Poly(polynomial + x, x, y), (x, y), (RationalSample(sp.S.Zero),), y
    )
    interval = RationalInterval(sp.S.One, sp.S(2))
    context._remember_root_rank_cert(1, interval, 0, 1)
    selected = AlgebraicRoot(sp.Poly(polynomial, y), interval, 1, 1, None, context)
    other = AlgebraicRoot(
        sp.Poly(2 * y - 3, y),
        RationalInterval(sp.Rational(3, 2), sp.Rational(3, 2)),
        0,
        1,
        None,
        FiberRootContext(sp.Poly(2 * y - 3 + x, x, y), (x, y), (RationalSample(sp.S.Zero),), y),
    )
    assert compare_samples(selected, other) == 1  # 2 > 3/2.
    assert context.roots_left_of(sp.Rational(3, 2)) == 1


def test_specialization_changes_multiplicity():
    from semialg.algebraic.roots import isolate_real_roots

    clear_algebraic_caches()
    x, y = sp.symbols("x y", real=True)
    alpha = isolate_real_roots(sp.Poly(x**2 - 2, x))[-1]
    polynomial = (y - 2) ** 2 * (y + 1)
    context = FiberRootContext(sp.Poly(polynomial + x**2 - 2, x, y), (x, y), (alpha,), y)
    interval = RationalInterval(sp.S.One, sp.S(3))
    context._remember_root_rank_cert(1, interval, 0, 1)
    selected = AlgebraicRoot(sp.Poly(polynomial, y), interval, 1, 1, None, context)
    other = AlgebraicRoot(
        sp.Poly(10 * y - 21, y),
        RationalInterval(sp.S(2), sp.S(4)),
        0,
        1,
        None,
        FiberRootContext(sp.Poly(10 * y - 21 + x**2 - 2, x, y), (x, y), (alpha,), y),
    )
    assert compare_samples(selected, other) == -1  # 2 < 21/10.
