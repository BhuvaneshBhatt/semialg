from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ..cad import PROJECTION_COLLINS, build_projection_sets, cad_decompose_full
from ..model import FallbackGuarantee, ProjectionConfig


@dataclass
class CollinsCompleteBackend:
    name: str = "collins_complete"

    def projection(
        self,
        polys: Sequence[sp.Expr],
        vars_: Sequence[sp.Symbol],
        config: ProjectionConfig | None = None,
    ):
        cfg = config or ProjectionConfig(
            operator=PROJECTION_COLLINS, strict_reduced_proj=False, strict_lift_conds=False
        )
        return build_projection_sets(polys, vars_, mode=PROJECTION_COLLINS, config=cfg)

    def decompose(
        self,
        polys: Sequence[sp.Expr],
        vars_: Sequence[sp.Symbol],
        config: ProjectionConfig | None = None,
        ecs_by_level: dict[int, Sequence[sp.Expr]] | None = None,
    ):
        cfg = config or ProjectionConfig(
            operator=PROJECTION_COLLINS, strict_reduced_proj=False, strict_lift_conds=False
        )
        cad = cad_decompose_full(
            polys,
            vars_,
            projection_mode=PROJECTION_COLLINS,
            projection_config=cfg,
            ecs_by_level=ecs_by_level or {},
        )
        cad.fallback_guarantee = FallbackGuarantee(
            globally_safe=True,
            reason="full Collins projection/lifting path",
            fallback_backend="collins_complete",
        )
        return cad
