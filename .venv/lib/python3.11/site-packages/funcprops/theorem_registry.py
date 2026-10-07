"""Authoritative declarative theorem registry for cross-property reasoning."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from .results import Boundedness, Continuity, Convexity, Monotonicity, TruthValue


@dataclass(frozen=True)
class TheoremRecord:
    theorem_id: str
    premise_property: str
    premise_value: object
    conclusions: tuple[tuple[str, object], ...]
    proof_basis: str
    completeness: str = "complete-on-same-domain"
    mutation_tags: tuple[str, ...] = ()
    verifier: Callable[[object], bool] | None = None


THEOREMS = (
    TheoremRecord(
        "calc.strict-increasing.injective",
        "monotonicity",
        Monotonicity.STRICTLY_INCREASING,
        (("injective", TruthValue.TRUE),),
        "strict monotonicity implies injectivity",
        mutation_tags=("weaken-strict", "change-domain"),
    ),
    TheoremRecord(
        "calc.strict-decreasing.injective",
        "monotonicity",
        Monotonicity.STRICTLY_DECREASING,
        (("injective", TruthValue.TRUE),),
        "strict monotonicity implies injectivity",
        mutation_tags=("weaken-strict", "change-domain"),
    ),
    TheoremRecord(
        "convex.strong.strict",
        "convexity",
        Convexity.STRONGLY_CONVEX,
        (("convexity", Convexity.STRICTLY_CONVEX), ("convexity", Convexity.CONVEX)),
        "strong convexity implication",
        mutation_tags=("drop-modulus",),
    ),
    TheoremRecord(
        "convex.strict.convex",
        "convexity",
        Convexity.STRICTLY_CONVEX,
        (("convexity", Convexity.CONVEX),),
        "strict convexity implication",
    ),
    TheoremRecord(
        "bounded.two-sided",
        "boundedness",
        Boundedness.BOUNDED,
        (("bounded_above", TruthValue.TRUE), ("bounded_below", TruthValue.TRUE)),
        "definition of boundedness",
    ),
    TheoremRecord(
        "continuous.truth",
        "continuity",
        Continuity.CONTINUOUS,
        (("continuous", TruthValue.TRUE),),
        "continuity classification projection",
    ),
)
_BY_PREMISE = {(t.premise_property, t.premise_value): t for t in THEOREMS}


def matching_theorem(prop, value):
    return _BY_PREMISE.get((prop, value))


def validate_registry():
    ids = [t.theorem_id for t in THEOREMS]
    return len(ids) == len(set(ids)) and all(
        t.conclusions and t.proof_basis and t.completeness for t in THEOREMS
    )


def generated_adversarial_cases():
    """Metadata-driven near misses used by architecture/adversarial tests."""
    return tuple((t.theorem_id, tag) for t in THEOREMS for tag in t.mutation_tags)
