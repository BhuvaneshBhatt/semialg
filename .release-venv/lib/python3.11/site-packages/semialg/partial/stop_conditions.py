from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class StopDecision:
    quantifier: str
    child_truth: bool
    should_stop: bool


def quantifier_stop_decision(quantifier: str, child_truth: bool) -> StopDecision:
    quantifier = quantifier.lower()
    if quantifier == "exists":
        return StopDecision(
            quantifier=quantifier, child_truth=child_truth, should_stop=bool(child_truth)
        )
    if quantifier == "forall":
        return StopDecision(
            quantifier=quantifier, child_truth=child_truth, should_stop=not bool(child_truth)
        )
    raise ValueError(f"Unsupported quantifier: {quantifier}")


__all__ = ["StopDecision", "quantifier_stop_decision"]
