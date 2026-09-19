from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import sympy as sp


@dataclass
class CanonIntSolveResult:
    variables: tuple[sp.Symbol, ...]
    formula: sp.Expr
    solutions: list[tuple[sp.Expr, ...]] = field(default_factory=list)
    method: str = "unknown"
    complete: bool = False
    provenance: list[str] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)


def dedup_int_points(points: Iterable[Sequence[sp.Expr]]) -> list[tuple[sp.Expr, ...]]:
    seen: set[tuple[str, ...]] = set()
    out: list[tuple[sp.Expr, ...]] = []
    for pt in points:
        t = tuple(sp.simplify(v) for v in pt)
        key = tuple(sp.srepr(v) for v in t)
        if key not in seen:
            seen.add(key)
            out.append(t)
    return out


def finite_points_to_formula(
    variables: Sequence[sp.Symbol], points: Iterable[Sequence[sp.Expr]]
) -> sp.Expr:
    pts = dedup_int_points(points)
    if not pts:
        return sp.false
    clauses = [sp.And(*[sp.Eq(v, val) for v, val in zip(variables, pt, strict=True)]) for pt in pts]
    return sp.Or(*clauses) if len(clauses) > 1 else clauses[0]


def canon_int_result(
    variables: Sequence[sp.Symbol],
    *,
    formula: sp.Expr | None = None,
    solutions: Iterable[Sequence[sp.Expr]] = (),
    method: str = "unknown",
    complete: bool = False,
    provenance: Iterable[str] = (),
    metadata: dict | None = None,
) -> CanonIntSolveResult:
    vars_t = tuple(variables)
    sols = dedup_int_points(solutions)
    form = sp.simplify(formula if formula is not None else finite_points_to_formula(vars_t, sols))
    return CanonIntSolveResult(
        variables=vars_t,
        formula=form,
        solutions=sols,
        method=method,
        complete=complete,
        provenance=list(provenance),
        metadata=dict(metadata or {}),
    )


__all__ = [
    "CanonIntSolveResult",
    "dedup_int_points",
    "finite_points_to_formula",
    "canon_int_result",
]
