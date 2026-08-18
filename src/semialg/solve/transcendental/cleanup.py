from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import sympy as sp


@dataclass(frozen=True)
class CleanupResult:
    original: sp.Expr
    cleaned: sp.Expr
    removed_redundancies: bool = False
    method: str = "transcendental_cleanup"
    metadata: dict = field(default_factory=dict)


def finite_points_form(
    variables: Sequence[sp.Symbol], points: Sequence[Sequence[sp.Expr]]
) -> sp.Expr:
    variables = tuple(variables)
    pts = sorted(set(tuple(p) for p in points), key=sp.default_sort_key)
    if not pts:
        return sp.false
    if len(variables) == 1:
        v = variables[0]
        values = [pt[0] for pt in pts]
        clauses = [sp.Eq(v, a) for a in values]
        return sp.Or(*clauses) if len(clauses) > 1 else clauses[0]
    clauses = [sp.And(*[sp.Eq(v, a) for v, a in zip(variables, pt, strict=True)]) for pt in pts]
    return sp.Or(*clauses) if len(clauses) > 1 else clauses[0]


def remove_redundant_disjunc(expr: sp.Expr) -> CleanupResult:
    if not isinstance(expr, sp.Or):
        cleaned = sp.simplify(expr)
        return CleanupResult(original=expr, cleaned=cleaned, removed_redundancies=(cleaned != expr))
    args = []
    for arg in expr.args:
        if arg not in args:
            args.append(arg)
    cleaned = sp.simplify(sp.Or(*args)) if len(args) > 1 else args[0]
    return CleanupResult(
        original=expr, cleaned=cleaned, removed_redundancies=(len(args) != len(expr.args))
    )


def recon_solved_points(
    variables: Sequence[sp.Symbol], points: Sequence[Sequence[sp.Expr]]
) -> CleanupResult:
    original = finite_points_form(variables, points)
    cleaned = sp.simplify(original)
    return CleanupResult(
        original=original,
        cleaned=cleaned,
        removed_redundancies=(cleaned != original),
        method="solved_form_reconstruction",
        metadata={"point_count": len(list(points))},
    )


def recon_univar_intv_form(
    variable: sp.Symbol, intervals: Sequence[tuple[sp.Expr, sp.Expr]]
) -> CleanupResult:
    clauses = []
    for a, b in intervals:
        clauses.append(sp.And(variable > a, variable < b))
    original = sp.Or(*clauses) if len(clauses) > 1 else (clauses[0] if clauses else sp.false)
    cleaned = sp.simplify(original)
    return CleanupResult(
        original=original,
        cleaned=cleaned,
        removed_redundancies=(cleaned != original),
        method="univariate_interval_reconstruction",
        metadata={"interval_count": len(tuple(intervals))},
    )


__all__ = [
    "CleanupResult",
    "finite_points_form",
    "remove_redundant_disjunc",
    "recon_solved_points",
    "recon_univar_intv_form",
]
