from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ..ec.propagation import propagate_ecs_by_family
from ..formula import Formula
from ..tti import extract_formula_families


@dataclass(frozen=True)
class TTIFamilySummary:
    index: int
    polynomial_count: int
    designated_ec: sp.Expr | None
    levels_with_ecs: tuple[int, ...]


def summarize_form_fams(
    formula: Formula, vars_: Sequence[sp.Symbol]
) -> tuple[TTIFamilySummary, ...]:
    families = extract_formula_families(formula)
    fam_maps = propagate_ecs_by_family(formula, vars_)
    out = []
    for idx, (fam, fam_map) in enumerate(zip(families, fam_maps, strict=True)):
        out.append(
            TTIFamilySummary(
                index=idx,
                polynomial_count=len(fam.polynomials),
                designated_ec=fam_map.designated,
                levels_with_ecs=tuple(sorted(fam_map.by_level)),
            )
        )
    return tuple(out)


__all__ = ["TTIFamilySummary", "summarize_form_fams"]
