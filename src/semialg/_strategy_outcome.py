"""Explicit outcomes for exact specialized-algorithm dispatch."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Generic, TypeVar

T = TypeVar("T")


class StrategyStatus(StrEnum):
    SUCCESS = "success"
    NOT_APPLICABLE = "not_applicable"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class StrategyOutcome(Generic[T]):
    status: StrategyStatus
    value: T | None = None
    reason: str | None = None

    @classmethod
    def success(cls, value: T):
        return cls(StrategyStatus.SUCCESS, value)

    @classmethod
    def not_applicable(cls, reason: str | None = None):
        return cls(StrategyStatus.NOT_APPLICABLE, reason=reason)

    @classmethod
    def unknown(cls, reason: str | None = None):
        return cls(StrategyStatus.UNKNOWN, reason=reason)
