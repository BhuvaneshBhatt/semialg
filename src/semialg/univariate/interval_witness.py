from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ..validation.solution_checking import form_sat_by_assign


@dataclass(frozen=True)
class IntervalWitnessCandidate:
    interval: tuple[object, object]
    witness: object | None


def witness_from_interval(interval: tuple[object, object], *, prefer_rational: bool = True):
    left, right = map(sp.sympify, interval)
    if left == right:
        return left
    midpoint = sp.simplify((left + right) / 2)
    return sp.nsimplify(midpoint) if prefer_rational else midpoint


def intv_wits_for_form(
    formula: sp.Expr,
    variable: sp.Symbol,
    intervals: Sequence[tuple[object, object]],
    *,
    domain=sp.Reals,
) -> list[IntervalWitnessCandidate]:
    out = []
    for interval in intervals:
        witness = witness_from_interval(interval)
        ok = form_sat_by_assign(
            formula, {variable: witness}, domain=domain, check_numeric_equalities=True
        )
        out.append(IntervalWitnessCandidate(interval, witness if ok else None))
    return out


__all__ = ["IntervalWitnessCandidate", "witness_from_interval", "intv_wits_for_form"]
