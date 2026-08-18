from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import sympy as sp

NumberLike = sp.Expr


@dataclass(frozen=True)
class Cell:
    level: int
    index: tuple[int, ...]
    sample: tuple[NumberLike, ...]
    intervals: tuple[tuple[NumberLike | None, NumberLike | None], ...]
    parent_index: tuple[int, ...] | None = None

    def describe(self, vars_: Sequence[sp.Symbol]) -> str:
        parts: list[str] = []
        for i, var in enumerate(vars_[: self.level]):
            left, right = self.intervals[i]
            if left is None and right is None:
                region = "(-oo, oo)"
            elif left is None:
                region = f"(-oo, {sp.sstr(right)})"
            elif right is None:
                region = f"({sp.sstr(left)}, oo)"
            elif sp.simplify(left - right) == 0:
                region = f"{{{sp.sstr(left)}}}"
            else:
                region = f"({sp.sstr(left)}, {sp.sstr(right)})"
            parts.append(f"{var}: {region}")
        return " | ".join(parts) if parts else "root"


@dataclass(frozen=True)
class ProjectionInfo:
    level: int
    variable: sp.Symbol | None
    primary_exprs: tuple[sp.Expr, ...] = field(default_factory=tuple)
    equational_constraints: tuple[sp.Expr, ...] = field(default_factory=tuple)
    mode: str = "collins"
    reduced_proj_used: bool = False
    designated_ec: sp.Expr | None = None
    family_ecs: tuple[sp.Expr, ...] = field(default_factory=tuple)
    lifting_exprs: tuple[sp.Expr, ...] = field(default_factory=tuple)
    strict_ok: bool = True
    lifting_strategy: str = "standard"
    tti_family_count: int = 0
    ec_selection_policy: str = "first"


@dataclass(frozen=True)
class NullificationEvent:
    level: int
    cell_index: tuple[int, ...]
    variable: sp.Symbol
    polynomial: sp.Expr
    delineating_polynomial: sp.Expr | None = None
    fallback_mode: str | None = None


@dataclass(frozen=True)
class DiagnosticEvent:
    kind: str
    message: str
    level: int | None = None
    cell_index: tuple[int, ...] | None = None
    details: dict[str, str] = field(default_factory=dict)


@dataclass
class CacheStats:
    hits: int = 0
    misses: int = 0


@dataclass
class Diagnostics:
    events: list[DiagnosticEvent] = field(default_factory=list)
    cache_stats: dict[str, CacheStats] = field(default_factory=dict)
    counters: dict[str, int] = field(default_factory=dict)

    def record(
        self,
        kind: str,
        message: str,
        *,
        level: int | None = None,
        cell_index: tuple[int, ...] | None = None,
        **details: object,
    ) -> None:
        self.events.append(
            DiagnosticEvent(
                kind=kind,
                message=message,
                level=level,
                cell_index=cell_index,
                details={key: str(value) for key, value in details.items()},
            )
        )
        self.counters[kind] = self.counters.get(kind, 0) + 1

    def cache_hit(self, cache_name: str) -> None:
        stats = self.cache_stats.setdefault(cache_name, CacheStats())
        stats.hits += 1

    def cache_miss(self, cache_name: str) -> None:
        stats = self.cache_stats.setdefault(cache_name, CacheStats())
        stats.misses += 1


@dataclass(frozen=True)
class LiftingCertificate:
    level: int
    parent_index: tuple[int, ...]
    variable: sp.Symbol
    strategy: str
    designated_ec: sp.Expr | None = None
    used_full_stack: bool = False
    used_collins_fallback: bool = False
    nullified_polynomials: tuple[sp.Expr, ...] = field(default_factory=tuple)
    delineating_polynomials: tuple[sp.Expr, ...] = field(default_factory=tuple)
    proof_notes: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class FallbackGuarantee:
    globally_safe: bool
    reason: str = ""
    fallback_backend: str | None = None


@dataclass
class CADResult:
    projection_sets: dict[int, list[sp.Poly]]
    cells_by_level: dict[int, list[Cell]]
    children_by_parent: dict[tuple[int, ...], list[Cell]]
    root_cell: Cell
    vars: tuple[sp.Symbol, ...]
    truth_by_cell: dict[tuple[int, ...], bool] = field(default_factory=dict)
    pruned_parents: tuple[tuple[int, ...], ...] = field(default_factory=tuple)
    lazy_cells_by_level: dict[int, list[Cell]] = field(default_factory=dict)
    formula_pruned_cells: tuple[tuple[int, ...], ...] = field(default_factory=tuple)
    projection_mode: str = "collins"
    projection_metadata: dict[int, ProjectionInfo] = field(default_factory=dict)
    nullification_events: tuple[NullificationEvent, ...] = field(default_factory=tuple)
    well_oriented: bool = True
    used_fallback_projection: bool = False
    lazard_evaluation_used: bool = False
    diagnostics: Diagnostics = field(default_factory=Diagnostics)
    lifting_validated: bool = True
    lifting_supplemented: bool = False
    tticad_used: bool = False
    lifting_certificates: tuple[LiftingCertificate, ...] = field(default_factory=tuple)
    fallback_guarantee: FallbackGuarantee = field(
        default_factory=lambda: FallbackGuarantee(globally_safe=False, reason="not assessed")
    )


@dataclass(frozen=True)
class ProjectionConfig:
    operator: str = "collins"
    use_ecs: bool = True
    enable_well_orientedness: bool = True
    fallback_to_collins: bool = True
    use_min_delin_polys: bool = True
    strict_reduced_proj: bool = True
    strict_lift_conds: bool = True
    use_lazard_evaluation: bool = False
    lazard_strict_mode: bool = False
    diagnostics: bool = True
    use_tticad_projection: bool = True
    validate_lifting: bool = True
    supplement_lifting_roots: bool = True
    ec_selection_policy: str = "lowest_degree"
    tticad_cross_family_mode: str = "designated_only"
    refine_on_demand: bool = True


@dataclass(frozen=True)
class QEConfig:
    partial: bool = False
    truth_invariant: bool = True
    projection: ProjectionConfig = field(default_factory=ProjectionConfig)
    simplify_output: bool = True
    auto_order: str | None = None


@dataclass
class QEResult:
    vars: tuple[sp.Symbol, ...]
    free_vars: tuple[sp.Symbol, ...]
    quantified_vars: tuple[sp.Symbol, ...]
    is_sentence: bool
    truth_value: bool | None
    free_level_cells: list[tuple[Cell, bool]] | None
    cad: CADResult
    qe_formula: sp.Expr | None = None
    cells_visited: int = 0
    partial_evaluation: bool = False
    guided_pruning: bool = False
