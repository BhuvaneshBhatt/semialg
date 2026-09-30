"""Uniform structured proof diagnostics for public semialg results."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class ProofDiagnostics:
    method: str
    steps: tuple[str, ...]
    used_complete_qe: bool
    used_critical_values: bool
    metadata: Mapping[str, object]


def _metadata_present(value: object) -> bool:
    """Return whether optional metadata carries information without truth-testing it."""
    if value is None or value is False or value is sp.false:
        return False
    if isinstance(value, str):
        return value.strip().lower() not in {"", "false", "none"}
    if isinstance(value, (tuple, list, dict, set, frozenset)):
        return bool(value)
    return True


def structured_proof_diagnostics(result) -> ProofDiagnostics:
    """Normalize heterogeneous result metadata into a stable proof trace."""
    method = str(getattr(result, "method", type(result).__name__))
    metadata = dict(getattr(result, "diagnostics", {}) or {})
    explicit = getattr(result, "proof_trace", ())
    steps = list(explicit if explicit else ())
    if not steps:
        if _metadata_present(metadata.get("back_substitutions")):
            steps.append("algebraic_back_substitution")
        if _metadata_present(metadata.get("algebraization")):
            steps.append(f"algebraization:{metadata['algebraization']}")
        if (
            "critical" in method
            or _metadata_present(metadata.get("critical_values"))
            or _metadata_present(metadata.get("critical_value_set"))
        ):
            steps.append("critical_value_geometry")
        if "optimization" in method:
            steps.append("optimization_bounds")
        if "qe" in method or "cad" in method:
            steps.append("quantifier_elimination")
        steps.append(method)
    used_qe = (
        any("qe" in x or "cad" in x or "projection" in x for x in steps)
        or "qe" in method
        or "cad" in method
    )
    used_critical = (
        any("critical" in step for step in steps)
        or "critical" in method
        or _metadata_present(metadata.get("critical_values"))
        or _metadata_present(metadata.get("critical_value_set"))
    )
    return ProofDiagnostics(method, tuple(dict.fromkeys(steps)), used_qe, used_critical, metadata)


__all__ = ["ProofDiagnostics", "structured_proof_diagnostics"]
