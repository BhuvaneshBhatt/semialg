"""Universal explanation/completeness projection."""

from __future__ import annotations

from dataclasses import dataclass

from .diagnostics import Completeness
from .results import as_detailed_result


@dataclass(frozen=True)
class ProofExplanation:
    conclusion: object
    theorem_ids: tuple[str, ...] = ()
    routes: tuple[str, ...] = ()
    established: tuple[str, ...] = ()
    unresolved: tuple[str, ...] = ()
    completeness: Completeness = Completeness.UNKNOWN
    certificate_status: str = "not-recorded"


def explain_proof(value, *, completeness=None):
    d = as_detailed_result(value)
    ids = tuple(e.theorem_id for e in d.evidence if getattr(e, "theorem_id", None))
    routes = (d.source,) if d.source else ()
    established = tuple(type(e).__name__ for e in d.evidence)
    unknown = d.primary is None or getattr(d.primary, "value", None) == "unknown"
    unresolved = ((d.source or "insufficient-evidence"),) if unknown else ()
    if completeness is None:
        if d.complete is True:
            comp = Completeness.COMPLETE_ON_DOMAIN
        elif d.complete is False:
            comp = Completeness.PARTIAL
        else:
            comp = Completeness.UNKNOWN
    else:
        comp = Completeness(completeness)
    cert = "not-recorded"
    if d.certificate is not None:
        v = getattr(d.certificate, "verify", None)
        cert = "verified" if callable(v) and v() else "present"
    return ProofExplanation(
        d.primary if d.primary is not None else "conditional/unknown",
        ids,
        routes,
        established,
        unresolved,
        comp,
        cert,
    )
