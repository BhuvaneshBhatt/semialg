"""Classification result types shared by all zero-testing stages."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Verdict(Enum):
    ZERO_PROVEN = "zero (proven)"
    ZERO_UNPROVEN = "zero-like (unproven)"
    NONZERO_PROVEN = "nonzero (proven)"
    NONZERO_LIKELY = "nonzero (likely)"
    UNKNOWN = "could not determine"


@dataclass(frozen=True)
class ZeroClassification:
    """A zero-test conclusion together with the evidence supporting it.

    ``is_zero`` is intentionally proof-oriented: it is ``True`` only for a
    proved zero, ``False`` for a proved or witnessed nonzero value, and
    ``None`` for zero-like numerical/probabilistic evidence or an unknown
    result.  Use ``suggests_zero`` when heuristic zero evidence is useful.
    """

    verdict: Verdict
    method: str
    detail: str = ""
    precision_bits: int | None = None
    trials: int | None = None
    evidence: str | None = None
    error_bound: float | None = None
    requested_error: float | None = None
    enclosure_history: tuple[str, ...] = ()

    @property
    def is_zero(self) -> bool | None:
        if self.verdict is Verdict.ZERO_PROVEN:
            return True
        if self.verdict in (Verdict.NONZERO_PROVEN, Verdict.NONZERO_LIKELY):
            return False
        return None

    @property
    def suggests_zero(self) -> bool:
        return self.verdict in (Verdict.ZERO_PROVEN, Verdict.ZERO_UNPROVEN)

    @property
    def proven(self) -> bool:
        return self.verdict in (Verdict.ZERO_PROVEN, Verdict.NONZERO_PROVEN)

    @property
    def certainty(self) -> str:
        """Return the strength of the evidence behind this classification.

        ``"certified"`` means the selected stage supplied an exact or
        rigorous certificate. ``"probable"`` denotes one-sided evidence
        such as a numerical witness. ``"heuristic"`` denotes suggestive
        evidence that is not strong enough for :func:`zerotest` to return a
        Boolean. ``"unknown"`` means no useful conclusion was reached.
        """
        if self.verdict in (Verdict.ZERO_PROVEN, Verdict.NONZERO_PROVEN):
            return "certified"
        if self.verdict is Verdict.NONZERO_LIKELY:
            return "probable"
        if self.verdict is Verdict.ZERO_UNPROVEN:
            return "heuristic"
        return "unknown"

    @property
    def reason(self) -> str:
        """Return a concise human-readable explanation of the decision."""
        parts = [self.method]
        if self.detail:
            parts.append(self.detail)
        if self.evidence:
            parts.append(f"evidence={self.evidence}")
        return ": ".join(parts)

    def __repr__(self):
        bits = f", {self.precision_bits} bits" if self.precision_bits else ""
        trials = f", {self.trials} trials" if self.trials else ""
        error = f", error<={self.error_bound:.3g}" if self.error_bound is not None else ""
        return (
            f"<ZeroClassification {self.verdict.value} via {self.method}"
            f"{bits}{trials}{error}: {self.detail}>"
        )
