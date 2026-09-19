"""Opt-in stage timing for the fast zero-testing oracle.

The normal :func:`zerotest` path does not construct a profiler.  Profiling is
therefore pay-for-play and intended for benchmarks, tuning, and diagnostics.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from time import perf_counter
from typing import TypeVar

from ._result import ZeroClassification

T = TypeVar("T")


@dataclass(frozen=True)
class StageTiming:
    """Timing and outcome metadata for one oracle stage."""

    stage: str
    seconds: float
    outcome: str
    detail: str = ""


@dataclass(frozen=True)
class ZeroTestProfile:
    """Result of an opt-in profiled zero-test run.

    The convenience properties mirror the final :class:`ZeroClassification`
    so callers can inspect why the oracle decided without unpacking internal
    timing records.
    """

    result: bool | None
    classification: ZeroClassification
    stages: tuple[StageTiming, ...]
    total_seconds: float

    @property
    def method(self) -> str:
        """Name of the stage that produced the final classification."""
        return self.classification.method

    @property
    def certainty(self) -> str:
        """Strength of the final evidence: certified, probable, heuristic, or unknown."""
        return self.classification.certainty

    @property
    def detail(self) -> str:
        """Stage-specific explanation of the final classification."""
        return self.classification.detail

    @property
    def evidence(self) -> str | None:
        """Optional compact description of the evidence used."""
        return self.classification.evidence

    @property
    def reason(self) -> str:
        """Concise human-readable explanation of why the oracle decided."""
        return self.classification.reason


class StageProfiler:
    """Small mutable recorder used only by ``profile_zerotest``."""

    __slots__ = ("_stages", "_start")

    def __init__(self) -> None:
        self._stages: list[StageTiming] = []
        self._start = perf_counter()

    def run(self, stage: str, func: Callable[[], T],
            outcome: Callable[[T], str] | None = None) -> T:
        start = perf_counter()
        value = func()
        elapsed = perf_counter() - start
        label = outcome(value) if outcome is not None else "done"
        self._stages.append(StageTiming(stage, elapsed, label))
        return value

    def note(self, stage: str, outcome: str, detail: str = "") -> None:
        self._stages.append(StageTiming(stage, 0.0, outcome, detail))

    def finish(
        self, result: bool | None, classification: ZeroClassification
    ) -> ZeroTestProfile:
        """Finalize timings together with the oracle's proof classification."""
        return ZeroTestProfile(
            result, classification, tuple(self._stages), perf_counter() - self._start
        )
