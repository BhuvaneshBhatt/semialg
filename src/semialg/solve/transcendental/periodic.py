from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp
from sympy.calculus.util import periodicity

from semialg.quantifiers import Exists

from .state import TransProblemState


@dataclass(frozen=True)
class PeriodicBoundingResult:
    variable: sp.Symbol
    period: sp.Expr | None
    lower_bound: sp.Expr
    upper_bound: sp.Expr
    assumptions: sp.Expr = sp.true


def detect_real_period(expr: sp.Expr, variable: sp.Symbol) -> sp.Expr | None:
    try:
        per = periodicity(expr, variable)
    except Exception:
        return None
    if per in (None, sp.S.ComplexInfinity):
        return None
    return sp.simplify(per)


def compute_periodic_window(expr: sp.Expr, variable: sp.Symbol) -> PeriodicBoundingResult | None:
    per = detect_real_period(expr, variable)
    if per is None:
        return None
    if per.is_real is False:
        return None
    return PeriodicBoundingResult(
        variable=variable,
        period=per,
        lower_bound=sp.Integer(0),
        upper_bound=sp.simplify(per),
        assumptions=sp.true,
    )


def _period_index(variable: sp.Symbol, *expressions: sp.Expr) -> sp.Symbol:
    """Create a non-colliding integer-index variable for periodic formulas."""

    sympified = tuple(sp.sympify(expr) for expr in expressions)
    occupied = set().union(*(expr.free_symbols for expr in sympified))
    occupied.add(variable)
    base = f"k_{variable.name}"
    index = sp.Symbol(base, real=True)
    counter = 1
    while index in occupied:
        index = sp.Symbol(f"{base}_{counter}", real=True)
        counter += 1
    return index


def periodic_intv_form(
    variable: sp.Symbol,
    representative_formula: sp.Expr,
    period: sp.Expr,
) -> sp.Expr:
    """Repeat a representative Boolean formula over all integer periods."""

    k = _period_index(variable, representative_formula, period)
    shifted = representative_formula.subs(variable, variable - k * period)
    return Exists(k, sp.And(sp.Contains(k, sp.S.Integers), shifted))


def recon_periodic_represent(
    variable: sp.Symbol,
    representatives: Sequence[sp.Expr],
    period: sp.Expr,
) -> sp.Expr:
    reps = tuple(sorted(set(sp.simplify(r) for r in representatives), key=sp.default_sort_key))
    if not reps:
        return sp.false
    k = _period_index(variable, period, *reps)
    root_formula = sp.Or(*(sp.Eq(variable, sp.simplify(root + k * period)) for root in reps))
    return Exists(k, sp.And(sp.Contains(k, sp.S.Integers), root_formula))


def recon_periodic_domain(
    variable: sp.Symbol,
    true_intervals: Sequence[tuple[sp.Expr, sp.Expr]],
    period: sp.Expr,
) -> sp.Expr:
    if not true_intervals:
        return sp.false
    endpoints = tuple(value for interval in true_intervals for value in interval)
    k = _period_index(variable, period, *endpoints)
    representative = sp.simplify(variable - k * period)
    clauses = [
        sp.And(representative > sp.simplify(a), representative < sp.simplify(b))
        for a, b in true_intervals
    ]
    interval_formula = sp.Or(*clauses) if len(clauses) > 1 else clauses[0]
    return Exists(k, sp.And(sp.Contains(k, sp.S.Integers), interval_formula))


def find_periodic_variables(state: TransProblemState) -> tuple[PeriodicBoundingResult, ...]:
    out = []
    for var in state.all_variables:
        box = compute_periodic_window(state.formula, var)
        if box is not None:
            out.append(box)
    return tuple(out)


__all__ = [
    "PeriodicBoundingResult",
    "detect_real_period",
    "compute_periodic_window",
    "periodic_intv_form",
    "recon_periodic_represent",
    "recon_periodic_domain",
    "find_periodic_variables",
]
