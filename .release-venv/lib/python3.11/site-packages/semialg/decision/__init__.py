"""Decision and structured semialgebraic solving APIs."""

from __future__ import annotations

from .api import (
    EquivalenceResult,
    ImplicationResult,
    IntervalComponent,
    SatisfiabilityResult,
    SemialgebraicSolution,
    TautologyResult,
    canonicalize_one_dimensional_formula,
    equivalent,
    implies,
    is_satisfiable,
    is_tautology,
    solve_semialgebraic,
)

__all__ = [
    "IntervalComponent",
    "SemialgebraicSolution",
    "EquivalenceResult",
    "ImplicationResult",
    "TautologyResult",
    "SatisfiabilityResult",
    "equivalent",
    "implies",
    "is_satisfiable",
    "is_tautology",
    "canonicalize_one_dimensional_formula",
    "solve_semialgebraic",
]
