"""Logical closure among already-certified function properties."""

from __future__ import annotations

from dataclasses import dataclass

from .theorem_registry import matching_theorem


@dataclass(frozen=True)
class TheoremFact:
    property: str
    value: object
    source: str = "given"
    theorem_id: str | None = None


def close_theorems(facts):
    out = {(f.property, f.value): f for f in facts}
    changed = True
    while changed:
        changed = False
        for key, _fact in tuple(out.items()):
            theorem = matching_theorem(*key)
            if theorem is None:
                continue
            for prop, value in theorem.conclusions:
                k = (prop, value)
                if k not in out:
                    out[k] = TheoremFact(
                        prop, value, f"theorem:{theorem.theorem_id}", theorem.theorem_id
                    )
                    changed = True
    return tuple(out.values())
