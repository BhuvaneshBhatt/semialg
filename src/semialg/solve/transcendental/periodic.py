from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp
from sympy.calculus.util import periodicity

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


def periodic_intv_form(
    variable: sp.Symbol,
    representative_formula: sp.Expr,
    period: sp.Expr,
) -> sp.Expr:
    k = sp.Symbol(f"k_{variable.name}", integer=True)
    lifted = sp.simplify(representative_formula.subs(variable, variable - k * period))
    return sp.Exists(k, lifted)


def recon_periodic_represent(
    variable: sp.Symbol,
    representatives: Sequence[sp.Expr],
    period: sp.Expr,
) -> sp.Expr:
    reps = tuple(sorted(set(sp.simplify(r) for r in representatives), key=sp.default_sort_key))
    if not reps:
        return sp.false
    k = sp.Symbol(f"k_{variable.name}", integer=True)
    clauses = [sp.Exists(k, sp.Eq(variable, sp.simplify(root + k * period))) for root in reps]
    return sp.Or(*clauses) if len(clauses) > 1 else clauses[0]


def recon_periodic_domain(
    variable: sp.Symbol,
    true_intervals: Sequence[tuple[sp.Expr, sp.Expr]],
    period: sp.Expr,
) -> sp.Expr:
    if not true_intervals:
        return sp.false
    k = sp.Symbol(f"k_{variable.name}", integer=True)
    clauses = []
    for a, b in true_intervals:
        clauses.append(
            sp.Exists(
                k,
                sp.And(
                    variable > sp.simplify(a + k * period), variable < sp.simplify(b + k * period)
                ),
            )
        )
    return sp.Or(*clauses) if len(clauses) > 1 else clauses[0]


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
