from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BackendSpec:
    name: str
    projection_family: str
    supports_ecs: bool
    supports_lazard: bool = False
    notes: str = ""


BACKENDS: dict[str, BackendSpec] = {
    "collins_complete": BackendSpec(
        "collins_complete", "full", False, False, "Certified-first Collins projection/lifting path."
    ),
    "collins": BackendSpec("collins", "full", False, False, "Classical conservative projection."),
    "mccallum": BackendSpec(
        "mccallum", "reduced", True, False, "McCallum-style reduced projection."
    ),
    "lazard": BackendSpec("lazard", "reduced", True, True, "Lazard-style lifting and evaluation."),
    "tticad": BackendSpec(
        "tticad", "tticad", True, False, "TTICAD-style family-aware reduced projection."
    ),
}


def get_backend(name: str) -> BackendSpec:
    try:
        return BACKENDS[name]
    except KeyError as exc:
        raise ValueError(f"Unsupported projection backend: {name}") from exc
