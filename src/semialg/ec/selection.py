from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ..formula import Formula, equational_constraints
from ..tti import extract_formula_families
from .scoring import ECScore, choose_eq_cons, rank_eq_cons


@dataclass(frozen=True)
class DesignatedECChoice:
    expr: sp.Expr | None
    policy: str
    ranked: tuple[ECScore, ...]


def choose_designated_ec(
    exprs: Sequence[sp.Expr], *, policy: str = "lowest_degree"
) -> DesignatedECChoice:
    ranked = rank_eq_cons(exprs)
    chosen = choose_eq_cons(exprs, policy=policy)
    return DesignatedECChoice(expr=chosen, policy=policy, ranked=ranked)


def choose_designated_form(
    formula: Formula, *, policy: str = "lowest_degree"
) -> DesignatedECChoice:
    return choose_designated_ec(equational_constraints(formula), policy=policy)


def choose_designated_fam(
    formula: Formula, *, policy: str = "lowest_degree"
) -> tuple[DesignatedECChoice, ...]:
    return tuple(
        choose_designated_ec(family.equational_constraints, policy=policy)
        for family in extract_formula_families(formula)
    )


__all__ = [
    "DesignatedECChoice",
    "choose_designated_ec",
    "choose_designated_form",
    "choose_designated_fam",
]
