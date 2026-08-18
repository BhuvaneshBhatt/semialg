from __future__ import annotations

from collections.abc import Iterable

from ..model import LiftingCertificate


def summarize_lifting_certs(certificates: Iterable[LiftingCertificate]) -> list[str]:
    lines: list[str] = []
    for cert in certificates:
        parts = [
            f"level={cert.level}",
            f"parent={cert.parent_index}",
            f"variable={cert.variable}",
            f"strategy={cert.strategy}",
        ]
        if cert.designated_ec is not None:
            parts.append("designated_ec=yes")
        if cert.used_collins_fallback:
            parts.append("collins_fallback=yes")
        if cert.used_full_stack:
            parts.append("full_stack=yes")
        if cert.nullified_polynomials:
            parts.append(f"nullifications={len(cert.nullified_polynomials)}")
        lines.append(", ".join(parts))
    return lines
