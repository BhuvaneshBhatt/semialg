from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp

from ..model import ProjectionInfo


@dataclass(frozen=True)
class ProjectionCertificate:
    level: int
    variable: sp.Symbol | None
    backend: str
    primary_exprs: tuple[sp.Expr, ...] = field(default_factory=tuple)
    designated_ec: sp.Expr | None = None
    strict_ok: bool = True
    notes: tuple[str, ...] = field(default_factory=tuple)


def build_proj_certs(metadata: dict[int, ProjectionInfo]) -> dict[int, ProjectionCertificate]:
    out: dict[int, ProjectionCertificate] = {}
    for level, info in metadata.items():
        notes = []
        if info.reduced_proj_used:
            notes.append("reduced projection applied")
        if info.strict_ok:
            notes.append("strict conditions satisfied")
        else:
            notes.append("strict conditions not fully certified")
        out[level] = ProjectionCertificate(
            level=level,
            variable=info.variable,
            backend=info.mode,
            primary_exprs=tuple(info.primary_exprs),
            designated_ec=info.designated_ec,
            strict_ok=info.strict_ok,
            notes=tuple(notes),
        )
    return out
