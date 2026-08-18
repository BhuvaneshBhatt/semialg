from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import cmp_to_key

import sympy as sp


@dataclass(frozen=True)
class Interval1D:
    """A one-dimensional interval with optional infinite endpoints.

    ``None`` denotes an unbounded endpoint. Finite endpoints are SymPy
    expressions, usually rational values or exact algebraic numbers emitted by
    the CAD sample layer.
    """

    left: sp.Expr | None
    right: sp.Expr | None
    left_closed: bool = False
    right_closed: bool = False

    @property
    def is_point(self) -> bool:
        return self.left is not None and self.right is not None and _equal(self.left, self.right)

    @property
    def is_universal(self) -> bool:
        return self.left is None and self.right is None


def _numeric(expr: sp.Expr | None, *, left: bool) -> sp.Expr:
    if expr is None:
        return -sp.oo if left else sp.oo
    return expr


def _compare_expr(a: sp.Expr | None, b: sp.Expr | None, *, left: bool = True) -> int:
    aa = _numeric(a, left=left)
    bb = _numeric(b, left=left)
    if aa is -sp.oo or aa == -sp.oo:
        return -1 if bb != -sp.oo else 0
    if bb is -sp.oo or bb == -sp.oo:
        return 1
    if aa is sp.oo or aa == sp.oo:
        return 1 if bb != sp.oo else 0
    if bb is sp.oo or bb == sp.oo:
        return -1
    diff = sp.simplify(aa - bb)
    if diff == 0:
        return 0
    try:
        sign = sp.sign(diff)
        if sign == -1:
            return -1
        if sign == 1:
            return 1
    except Exception:
        pass
    aval = sp.N(aa, 80)
    bval = sp.N(bb, 80)
    if aval < bval:
        return -1
    if aval > bval:
        return 1
    return 0


def _equal(a: sp.Expr, b: sp.Expr) -> bool:
    return _compare_expr(a, b) == 0


def _left_sort_cmp(a: Interval1D, b: Interval1D) -> int:
    cmp = _compare_expr(a.left, b.left, left=True)
    if cmp != 0:
        return cmp
    # closed left endpoints come first so [a,b] can absorb (a,c).
    if a.left_closed != b.left_closed:
        return -1 if a.left_closed else 1
    return _compare_expr(a.right, b.right, left=False)


def intvs_overlap_or_touch(left: Interval1D, right: Interval1D) -> bool:
    if left.right is None or right.left is None:
        return True
    cmp = _compare_expr(right.left, left.right)
    if cmp < 0:
        return True
    if cmp > 0:
        return False
    return left.right_closed or right.left_closed


def _right_max(left: Interval1D, right: Interval1D) -> tuple[sp.Expr | None, bool]:
    if left.right is None or right.right is None:
        return None, False
    cmp = _compare_expr(left.right, right.right, left=False)
    if cmp < 0:
        return right.right, right.right_closed
    if cmp > 0:
        return left.right, left.right_closed
    return left.right, left.right_closed or right.right_closed


def merge_intervals(intervals: Sequence[Interval1D]) -> tuple[Interval1D, ...]:
    """Merge overlapping or endpoint-touching intervals exactly where possible."""

    if not intervals:
        return tuple()
    ordered = sorted(intervals, key=cmp_to_key(_left_sort_cmp))
    merged: list[Interval1D] = []
    for interval in ordered:
        if not merged:
            merged.append(interval)
            continue
        prev = merged[-1]
        if not intvs_overlap_or_touch(prev, interval):
            merged.append(interval)
            continue
        right, right_closed = _right_max(prev, interval)
        merged[-1] = Interval1D(prev.left, right, prev.left_closed, right_closed)
    return tuple(merged)


def interval_condition(
    var: sp.Symbol,
    left: sp.Expr | None,
    right: sp.Expr | None,
    *,
    left_closed: bool = False,
    right_closed: bool = False,
) -> sp.Expr:
    pieces: list[sp.Expr] = []
    if left is not None and right is not None and _equal(left, right):
        if left_closed and right_closed:
            return sp.Eq(var, left)
        return sp.false
    if left is not None:
        pieces.append(var >= left if left_closed else var > left)
    if right is not None:
        pieces.append(var <= right if right_closed else var < right)
    if not pieces:
        return sp.true
    return sp.And(*pieces)


def intervals_to_formula(var: sp.Symbol, intervals: Sequence[Interval1D]) -> sp.Expr:
    merged = merge_intervals(intervals)
    if not merged:
        return sp.false
    if len(merged) == 1 and merged[0].is_universal:
        return sp.true
    expr = sp.Or(
        *[
            interval_condition(
                var,
                interval.left,
                interval.right,
                left_closed=interval.left_closed,
                right_closed=interval.right_closed,
            )
            for interval in merged
        ]
    )
    try:
        return sp.simplify_logic(expr, form="dnf")
    except Exception:
        return expr


__all__ = ["Interval1D", "interval_condition", "merge_intervals", "intervals_to_formula"]
