from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from ..model import CADResult


@dataclass(frozen=True)
class TTILiftingPlan:
    level: int
    variable: sp.Symbol | None
    designated_ecs: tuple[sp.Expr, ...]
    lifting_exprs: tuple[sp.Expr, ...]
    reduced: bool


def build_tticad_plan(cad: CADResult) -> tuple[TTILiftingPlan, ...]:
    plans = []
    for level in sorted(cad.projection_metadata):
        info = cad.projection_metadata[level]
        plans.append(
            TTILiftingPlan(
                level=level,
                variable=info.variable,
                designated_ecs=tuple(info.family_ecs or info.equational_constraints),
                lifting_exprs=tuple(info.lifting_exprs),
                reduced=bool(info.reduced_proj_used),
            )
        )
    return tuple(plans)


__all__ = ["TTILiftingPlan", "build_tticad_plan"]
