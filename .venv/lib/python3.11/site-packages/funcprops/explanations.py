"""Structured proof and UNKNOWN explanations."""

from __future__ import annotations

from dataclasses import dataclass

from .results import as_detailed_result


@dataclass(frozen=True)
class ResultExplanation:
    conclusion: object
    route: str | None
    established: tuple[str, ...] = ()
    unresolved: tuple[str, ...] = ()
    certificate: str | None = None

    def __str__(self):
        parts = [f"Result: {self.conclusion}"]
        if self.route:
            parts.append(f"Proof route: {self.route}")
        if self.established:
            parts.append("Established: " + "; ".join(self.established))
        if self.unresolved:
            parts.append("Unresolved: " + "; ".join(self.unresolved))
        if self.certificate:
            parts.append(f"Certificate: {self.certificate}")
        return "\n".join(parts)


def explain(value):
    d = as_detailed_result(value)
    established = []
    unresolved = []
    if d.evidence:
        established.append(f"{len(d.evidence)} evidence item(s)")
    if d.primary is None or getattr(d.primary, "value", None) == "unknown":
        unresolved.append((d.source or "insufficient evidence").replace("-", " "))
    cert = type(d.certificate).__name__ if d.certificate is not None else None
    return ResultExplanation(
        d.primary if d.primary is not None else "conditional/unknown",
        d.source,
        tuple(established),
        tuple(unresolved),
        cert,
    )
