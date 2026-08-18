from __future__ import annotations

from dataclasses import dataclass, field

from .features import ProblemFeatures


@dataclass(frozen=True)
class StrategyMemoryEntry:
    signature: str
    backend: str
    variable_order: tuple[str, ...]
    partial: bool
    success: bool
    runtime: float | None = None
    notes: str = ""


@dataclass
class StrategyMemory:
    entries: list[StrategyMemoryEntry] = field(default_factory=list)

    def add(self, entry: StrategyMemoryEntry) -> None:
        self.entries.append(entry)

    def best_backend_signature(self, signature: str) -> str | None:
        candidates = [e for e in self.entries if e.signature == signature and e.success]
        if not candidates:
            return None
        candidates.sort(key=lambda e: (float("inf") if e.runtime is None else e.runtime, e.backend))
        return candidates[0].backend

    def summary(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for entry in self.entries:
            out[entry.backend] = out.get(entry.backend, 0) + 1
        return out


def feature_signature(features: ProblemFeatures) -> str:
    return "|".join(
        [
            f"vars={features.variable_count}",
            f"deg={features.max_total_degree}",
            f"atoms={features.num_atoms}",
            f"alts={features.quantifier_alternations}",
            f"or={int(features.has_disjunction)}",
            f"ec={int(features.has_ecs)}",
        ]
    )
