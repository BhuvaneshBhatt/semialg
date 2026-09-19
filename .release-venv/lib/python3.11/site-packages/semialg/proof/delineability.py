from __future__ import annotations

from collections.abc import Iterable

from ..model import NullificationEvent


def likely_delineable(events: Iterable[NullificationEvent]) -> bool:
    return not any(True for _ in events)
