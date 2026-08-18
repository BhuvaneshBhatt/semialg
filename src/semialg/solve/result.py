from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .domains import SolveDomain


@dataclass
class SolveResult:
    method: str
    domain: SolveDomain
    result: Any
    normalized_text: str | None = None
    preprocess_changed: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)


__all__ = ["SolveResult"]
