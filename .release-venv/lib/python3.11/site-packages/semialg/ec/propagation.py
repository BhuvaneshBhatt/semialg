from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ..formula import Formula
from ..tti import extract_formula_families
from .selection import choose_designated_ec


@dataclass(frozen=True)
class FamilyECMap:
    family_index: int
    equational_constraints: tuple[sp.Expr, ...]
    designated: sp.Expr | None
    by_level: dict[int, tuple[sp.Expr, ...]]


def _ecs_by_level(
    exprs: Sequence[sp.Expr], vars_: Sequence[sp.Symbol]
) -> dict[int, tuple[sp.Expr, ...]]:
    by_level: dict[int, tuple[sp.Expr, ...]] = {}
    exprs = tuple(sp.expand(e) for e in exprs)
    for level in range(1, len(vars_) + 1):
        allowed = set(vars_[:level])
        current = tuple(e for e in exprs if e.free_symbols and e.free_symbols <= allowed)
        if current:
            by_level[level] = current
    return by_level


def propagate_ecs_by_family(
    formula: Formula, vars_: Sequence[sp.Symbol], *, policy: str = "lowest_degree"
) -> tuple[FamilyECMap, ...]:
    out = []
    for idx, family in enumerate(extract_formula_families(formula)):
        choice = choose_designated_ec(family.equational_constraints, policy=policy)
        out.append(
            FamilyECMap(
                family_index=idx,
                equational_constraints=tuple(family.equational_constraints),
                designated=choice.expr,
                by_level=_ecs_by_level(family.equational_constraints, vars_),
            )
        )
    return tuple(out)


def merged_fam_ecs_by_level(
    formula: Formula, vars_: Sequence[sp.Symbol], *, policy: str = "lowest_degree"
) -> dict[int, tuple[sp.Expr, ...]]:
    merged: dict[int, list[sp.Expr]] = {}
    for fam in propagate_ecs_by_family(formula, vars_, policy=policy):
        for level, exprs in fam.by_level.items():
            merged.setdefault(level, [])
            for expr in exprs:
                if expr not in merged[level]:
                    merged[level].append(expr)
    return {k: tuple(v) for k, v in merged.items()}


__all__ = ["FamilyECMap", "propagate_ecs_by_family", "merged_fam_ecs_by_level"]
