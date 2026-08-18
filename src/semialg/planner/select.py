from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp

from ..cad.constants import (
    PROJECTION_COLLINS,
    PROJECTION_LAZARD,
    PROJECTION_MCCALLUM,
    PROJECTION_TTICAD,
)
from ..formula import Formula, ParsedPrenexFormula, formula_polynomials
from ..model import ProjectionConfig, QEConfig
from .analyze import ProblemAnalysis, analyze_formula, analyze_parsed_formula
from .heuristics import choose_best_var_order
from .strategy_memory import StrategyMemory, feature_signature

PROJECTION_EC = PROJECTION_MCCALLUM


@dataclass(frozen=True)
class StrategySelection:
    backend: str
    variable_order: tuple[sp.Symbol, ...]
    partial: bool
    projection: ProjectionConfig
    notes: tuple[str, ...] = field(default_factory=tuple)


def select_strat_analysis(
    analysis: ProblemAnalysis,
    *,
    strategy_memory: StrategyMemory | None = None,
) -> StrategySelection:
    f = analysis.features
    notes = list(analysis.notes)
    signature = feature_signature(f)
    memory_backend = (
        strategy_memory.best_backend_signature(signature) if strategy_memory is not None else None
    )

    backend = memory_backend or PROJECTION_COLLINS
    partial = False

    if memory_backend is not None:
        notes.append(f"Using memory-preferred backend: {memory_backend}.")
    elif f.is_univariate:
        backend = PROJECTION_COLLINS
        partial = False
        notes.append("Using conservative complete backend for univariate problem.")
    elif True:
        if f.has_disjunction:
            backend = PROJECTION_TTICAD
            partial = True
            notes.append("Selected TTICAD due to disjunction/branch structure.")
        elif f.has_ecs and f.variable_count >= 2:
            backend = PROJECTION_EC if f.ec_density >= 0.3 else PROJECTION_MCCALLUM
            partial = f.quantifier_alternations > 0
            notes.append(
                "Selected reduced projection because equational constraints are available."
            )
        elif f.quantifier_alternations > 0:
            backend = PROJECTION_MCCALLUM
            partial = True
            notes.append("Selected partial reduced CAD for alternating quantifiers.")
        elif f.max_total_degree >= 5:
            backend = PROJECTION_LAZARD
            partial = True
            notes.append("Selected Lazard-style mode for higher-degree algebraic behavior.")
        else:
            backend = PROJECTION_COLLINS
            partial = False
            notes.append("Selected conservative complete backend.")

    projection = ProjectionConfig(
        operator=backend,
        use_ecs=f.has_ecs,
        use_tticad_projection=(backend == PROJECTION_TTICAD),
    )
    return StrategySelection(
        backend=backend,
        variable_order=f.suggested_variable_order or f.variables,
        partial=partial,
        projection=projection,
        notes=tuple(notes),
    )


def select_strat_for_form(
    formula: Formula,
    *,
    parsed: ParsedPrenexFormula | None = None,
    strategy_memory: StrategyMemory | None = None,
) -> StrategySelection:
    analysis = analyze_parsed_formula(parsed) if parsed is not None else analyze_formula(formula)
    polys = tuple(formula_polynomials(formula))
    order = (
        choose_best_var_order(analysis.features, polys) if polys else analysis.features.variables
    )
    base = select_strat_analysis(analysis, strategy_memory=strategy_memory)
    return StrategySelection(
        backend=base.backend,
        variable_order=order,
        partial=base.partial,
        projection=base.projection,
        notes=base.notes,
    )


def build_qe_config(selection: StrategySelection, *, simplify_output: bool = True) -> QEConfig:
    return QEConfig(
        partial=selection.partial,
        truth_invariant=(selection.backend == PROJECTION_TTICAD or selection.partial),
        projection=selection.projection,
        simplify_output=simplify_output,
        auto_order="planner",
    )
