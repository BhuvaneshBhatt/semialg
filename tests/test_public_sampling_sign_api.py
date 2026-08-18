from __future__ import annotations

import sympy as sp

from semialg import sample_point, sample_points, sign_at, sign_vector


def test_sign_at_rational_mapping_and_sequence() -> None:
    x, y = sp.symbols("x y")
    assert sign_at(x**2 + y - 1, {x: 0, y: 1}) == 0
    assert sign_at(x + y, {x: sp.Rational(1, 3), y: sp.Rational(1, 2)}) == 1
    assert sign_at(x - y, [sp.Rational(1, 3), sp.Rational(1, 2)], variables=[x, y]) == -1


def test_sign_vector_preserves_polynomial_order() -> None:
    x, y = sp.symbols("x y")
    signs = sign_vector([x, y, x + y, x**2 + y**2 - 1], {x: 0, y: 1})
    assert signs == (0, 1, 1, 0)


def test_sign_at_handles_algebraic_values_exactly() -> None:
    x = sp.symbols("x")
    root_two = sp.sqrt(2)
    assert sign_at(x**2 - 2, {x: root_two}) == 0
    assert sign_at(x - 1, {x: root_two}) == 1
    assert sign_at(x - 2, {x: root_two}) == -1


def test_sample_point_returns_satisfying_witness_for_disk() -> None:
    x, y = sp.symbols("x y", real=True)
    witness = sample_point(x**2 + y**2 < 1, [x, y])
    assert witness is not None
    assert bool((x**2 + y**2 < 1).subs(witness)) is True


def test_sample_point_returns_none_for_empty_region() -> None:
    x, y = sp.symbols("x y", real=True)
    assert sample_point(x**2 + y**2 < 0, [x, y]) is None


def test_sample_points_finds_requested_distinct_witnesses_when_available() -> None:
    x = sp.symbols("x", real=True)
    witnesses = sample_points(x**2 > 1, [x], count=2)
    assert len(witnesses) == 2
    assert all(bool((x**2 > 1).subs(witness)) for witness in witnesses)
    assert len({sp.simplify(witness[x]) for witness in witnesses}) == 2
