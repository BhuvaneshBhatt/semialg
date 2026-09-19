from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..formula import Formula, ParsedPrenexFormula
from .features import ProblemFeatures, extract_problem_features


@dataclass(frozen=True)
class ProblemAnalysis:
    features: ProblemFeatures
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def estimated_difficulty(self) -> str:
        f = self.features
        if f.variable_count <= 1:
            return "low"
        score = 0
        score += max(f.variable_count - 1, 0)
        score += f.quantifier_alternations * 2
        score += int(f.has_disjunction)
        score += int(f.has_negation)
        score += int(f.max_total_degree >= 4)
        score += int(f.num_atoms >= 6)
        if score <= 1:
            return "low"
        if score <= 4:
            return "medium"
        return "high"


def analyze_formula(
    formula: Formula,
    *,
    variables: Iterable[sp.Symbol] | None = None,
    quantifiers: Sequence[tuple[str, sp.Symbol]] = (),
) -> ProblemAnalysis:
    features = extract_problem_features(formula, variables=variables, quantifiers=quantifiers)
    notes: list[str] = []
    if features.is_univariate:
        notes.append("Univariate or effectively univariate problem; direct solving may beat CAD.")
    if features.has_ecs:
        notes.append("Formula contains equational constraints; reduced projection may help.")
    if features.has_disjunction:
        notes.append(
            "Disjunction detected; TTICAD or branch-local processing is likely beneficial."
        )
    if features.quantifier_alternations:
        notes.append("Alternating quantifiers detected; prefer partial/lazy CAD where possible.")
    if features.max_total_degree >= 4:
        notes.append("High total degree detected; algebraic root isolation may dominate run time.")
    if not notes:
        notes.append("No special structural warnings detected.")
    return ProblemAnalysis(features=features, notes=tuple(notes))


def analyze_parsed_formula(parsed: ParsedPrenexFormula) -> ProblemAnalysis:
    return analyze_formula(parsed.matrix, variables=parsed.vars, quantifiers=parsed.quantifiers)
