from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ..cad import PROJECTION_TTICAD, build_projection_sets
from ..formula import Formula, formula_polynomials
from ..model import ProjectionConfig, ProjectionInfo
from ..tti import extract_formula_families


@dataclass(frozen=True)
class TTIProjectionResult:
    projection_sets: dict[int, list[sp.Poly]]
    metadata: dict[int, ProjectionInfo]
    family_count: int


def build_tticad_projection(
    formula: Formula, vars_: Sequence[sp.Symbol], *, config: ProjectionConfig | None = None
) -> TTIProjectionResult:
    config = config or ProjectionConfig(operator=PROJECTION_TTICAD)
    families = extract_formula_families(formula)
    proj, meta = build_projection_sets(
        formula_polynomials(formula),
        tuple(vars_),
        mode=PROJECTION_TTICAD,
        ecs_by_level={},
        config=config,
        formula_families=families,
    )
    return TTIProjectionResult(projection_sets=proj, metadata=meta, family_count=len(families))


__all__ = ["TTIProjectionResult", "build_tticad_projection"]
