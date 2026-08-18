from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from .symmetric_difference import find_grid_witness


@dataclass(frozen=True)
class EquivCounterex:
    point: Mapping[str, object]
    left: bool
    right: bool


@dataclass(frozen=True)
class EquivalenceReport:
    equivalent: bool
    checked_points: int
    mismatches: tuple[EquivCounterex, ...] = field(default_factory=tuple)
    method: str = "grid"
    proof_attempted: bool = False


def sym_diff_empty(
    left: sp.Expr, right: sp.Expr, variables: Sequence[sp.Symbol]
) -> EquivalenceReport:
    check = find_grid_witness(left, right, variables)
    if check.equivalent_on_grid:
        return EquivalenceReport(equivalent=True, checked_points=check.checked_points)
    assert check.witness is not None
    mismatch = EquivCounterex(
        point={sp.sstr(sym): val for sym, val in check.witness.assignment.items()},
        left=check.witness.left_value,
        right=check.witness.right_value,
    )
    return EquivalenceReport(
        equivalent=False, checked_points=check.checked_points, mismatches=(mismatch,)
    )


__all__ = ["EquivCounterex", "EquivalenceReport", "sym_diff_empty"]
