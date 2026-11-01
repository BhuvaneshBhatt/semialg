from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp


@dataclass(frozen=True)
class EquivCounterex:
    point: Mapping[sp.Symbol, object]
    left: bool
    right: bool


@dataclass(frozen=True)
class EquivalenceReport:
    equivalent: bool
    checked_points: int
    mismatches: tuple[EquivCounterex, ...] = field(default_factory=tuple)
    method: str = "symmetric_difference"
    proof_attempted: bool = False


def sym_diff_empty(
    left: sp.Expr, right: sp.Expr, variables: Sequence[sp.Symbol]
) -> EquivalenceReport:
    """Decide exact real set equivalence, retaining the historical report type.

    Finite sampling is available separately as ``find_grid_witness``; absence
    of a grid witness is never evidence of exact equivalence here.
    """
    from ..decision import equivalent, is_satisfiable
    from ..simplify.atoms import normalize_atoms

    result = equivalent(
        normalize_atoms(sp.sympify(left)),
        normalize_atoms(sp.sympify(right)),
        variables,
        return_result=True,
    )
    mismatches = ()
    if result.counterexample is not None:
        point = dict(result.counterexample)
        left_value = is_satisfiable(sp.sympify(left).subs(point), ())
        right_value = is_satisfiable(sp.sympify(right).subs(point), ())
        mismatches = (EquivCounterex(point, bool(left_value), bool(right_value)),)
    return EquivalenceReport(
        equivalent=result.equivalent,
        checked_points=len(mismatches),
        mismatches=mismatches,
        method=result.method,
        proof_attempted=True,
    )


__all__ = ["EquivCounterex", "EquivalenceReport", "sym_diff_empty"]
