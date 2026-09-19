"""Reusable exact semantic comparisons for differential tests."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import sympy as sp

from semialg import is_equal


def assert_formulas_equivalent(
    left: sp.Expr,
    right: sp.Expr,
    variables: Sequence[sp.Symbol],
) -> None:
    """Require two formulas to define exactly the same real set."""

    assert is_equal(left, right, tuple(variables)), (
        f"formulas are not semantically equivalent:\nleft={left}\nright={right}"
    )


def canonical_points(
    points: Iterable[Sequence[object] | Mapping[sp.Symbol, object]],
    variables: Sequence[sp.Symbol],
) -> frozenset[tuple[sp.Expr, ...]]:
    """Normalize exact points independently of container/result representation."""

    variables = tuple(variables)
    normalized: set[tuple[sp.Expr, ...]] = set()
    for point in points:
        if isinstance(point, Mapping):
            coordinates = tuple(sp.simplify(point[var]) for var in variables)
        else:
            coordinates = tuple(sp.simplify(value) for value in point)
        normalized.add(coordinates)
    return frozenset(normalized)


def assert_exact_point_sets_equal(
    left: Iterable[Sequence[object] | Mapping[sp.Symbol, object]],
    right: Iterable[Sequence[object] | Mapping[sp.Symbol, object]],
    variables: Sequence[sp.Symbol],
) -> None:
    """Require two finite exact solution sets to contain the same points."""

    assert canonical_points(left, variables) == canonical_points(right, variables)
