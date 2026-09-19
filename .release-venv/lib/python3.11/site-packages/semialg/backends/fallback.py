from __future__ import annotations

from dataclasses import dataclass

from ..model import FallbackGuarantee


@dataclass(frozen=True)
class FallbackDecision:
    use_backend: str
    guarantee: FallbackGuarantee
    reason: str = ""


def collins_safe_fallback(
    reason: str = "reduced backend assumptions not certified",
) -> FallbackDecision:
    return FallbackDecision(
        use_backend="collins_complete",
        guarantee=FallbackGuarantee(
            globally_safe=True, reason=reason, fallback_backend="collins_complete"
        ),
        reason=reason,
    )


def maybe_fallback(strict_ok: bool, requested_backend: str) -> FallbackDecision | None:
    if strict_ok or requested_backend in {"collins", "collins_complete"}:
        return None
    return collins_safe_fallback(f"strict conditions failed for backend {requested_backend}")
