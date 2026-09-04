from __future__ import annotations

from .atoms import canonicalize_relation, normalize_atoms
from .boolean import simplify_boolean
from .cell_union import cell_to_interval_1d, cell_union_to_formula, cell_union_to_intervals
from .equality import simplify_equalities
from .formula import simplify_qe_formula, simplify_semialgebraic_formula
from .implication import (
    ImplicationMinimizationStats,
    clear_implication_minimization_stats,
    implication_minimization_stats,
)
from .intervals import Interval1D, interval_condition, intervals_to_formula, merge_intervals

canonicalize_qe_formula = simplify_qe_formula

__all__ = [
    "ImplicationMinimizationStats",
    "Interval1D",
    "canonicalize_qe_formula",
    "clear_implication_minimization_stats",
    "cell_to_interval_1d",
    "cell_union_to_formula",
    "cell_union_to_intervals",
    "implication_minimization_stats",
    "interval_condition",
    "intervals_to_formula",
    "merge_intervals",
    "canonicalize_relation",
    "normalize_atoms",
    "simplify_boolean",
    "simplify_equalities",
    "simplify_qe_formula",
    "simplify_semialgebraic_formula",
]
