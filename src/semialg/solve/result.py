from __future__ import annotations

from dataclasses import dataclass, field

from .domains import SolveDomain


@dataclass
class SolveResult:
    method: str
    domain: SolveDomain
    result: object
    normalized_text: str | None = None
    preprocess_changed: bool = False
    metadata: dict[str, object] = field(default_factory=dict)


__all__ = ["SolveResult"]
