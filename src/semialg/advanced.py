from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class ProjectionTheory(str, Enum):
    COLLINS = "collins"
    MCCALLUM = "mccallum"
    BROWN = "brown"
    LAZARD = "lazard"


@dataclass(frozen=True)
class CompletenessStatus:
    theory: ProjectionTheory
    correctness_complete: bool
    reason: str


class UnsupCompleteError(NotImplementedError):
    """Raised when a user requests a formally correctness-complete backend that is not implemented."""


def get_completeness_status(theory: str | ProjectionTheory) -> CompletenessStatus:
    theory = ProjectionTheory(theory)
    if theory == ProjectionTheory.COLLINS:
        return CompletenessStatus(
            theory=theory,
            correctness_complete=False,
            reason=(
                "The current package contains a Collins-style implementation, "
                "not a formally verified or research-grade complete implementation."
            ),
        )
    if theory == ProjectionTheory.MCCALLUM:
        return CompletenessStatus(
            theory=theory,
            correctness_complete=False,
            reason=(
                "A full McCallum implementation requires complete well-orientedness theory, nullification "
                "recovery, reduced lifting, and correctness conditions that are not fully implemented here."
            ),
        )
    if theory == ProjectionTheory.BROWN:
        return CompletenessStatus(
            theory=theory,
            correctness_complete=False,
            reason=(
                "A full Brown implementation requires the exact Brown reduced projection and its lifting "
                "conditions, which are not fully implemented in this package."
            ),
        )
    if theory == ProjectionTheory.LAZARD:
        return CompletenessStatus(
            theory=theory,
            correctness_complete=False,
            reason=(
                "A full Lazard implementation requires Lazard valuation, Lazard evaluation during lifting, "
                "and the full validated correctness machinery; this package does not yet provide that."
            ),
        )
    raise AssertionError("unreachable")


def require_correct_complete(theory: str | ProjectionTheory) -> None:
    status = get_completeness_status(theory)
    if not status.correctness_complete:
        raise UnsupCompleteError(status.reason)


def explain_missing_comps() -> dict[str, list[str]]:
    return {
        "mccallum": [
            "Formal reduced projection operator rather than the current approximation.",
            "Complete well-orientedness tests and nullification handling.",
            "Correct reduced lifting rules and equational-constraint lifting restrictions.",
            "Exhaustive regression corpus against trusted CAD systems.",
        ],
        "brown": [
            "The exact Brown projection refinements, not just McCallum-style approximations.",
            "Brown-specific correctness conditions for lifting.",
            "Validated interaction with equational constraints and pruning.",
        ],
        "lazard": [
            "Lazard valuation and valuation-invariant lifting.",
            "Lazard evaluation at sample points with repeated factor cancellation.",
            "The full validity-proof-aligned implementation, including tricky degenerate cases.",
        ],
    }
