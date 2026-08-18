from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

from .cases import BenchmarkCase


@dataclass(frozen=True)
class RunMetrics:
    elapsed_seconds: float = 0.0
    cells_total: int = 0
    top_level_cells: int = 0
    cells_visited: int = 0
    diagnostics_counters: Mapping[str, int] = field(default_factory=dict)


@dataclass(frozen=True)
class BenchmarkRunSpec:
    name: str
    config: object | None = None


@dataclass
class BenchmarkRunResult:
    spec: BenchmarkRunSpec
    qe_result: object | None = None
    metrics: RunMetrics = field(default_factory=RunMetrics)
    sample_truth: tuple[bool, ...] = ()


@dataclass(frozen=True)
class DifferentialMismatch:
    point: Mapping[str, object]
    left: bool
    right: bool
    left_run: str
    right_run: str


@dataclass
class BenchmarkCaseResult:
    case: BenchmarkCase
    parsed: object
    runs: dict[str, BenchmarkRunResult]
    checker_mismatches: list[tuple[str, Mapping[str, object], bool, bool]] = field(
        default_factory=list
    )
    differential_mismatches: list[DifferentialMismatch] = field(default_factory=list)
    equivalence_reports: dict[str, object] = field(default_factory=dict)
    taxonomy: tuple[str, ...] = field(default_factory=tuple)

    @property
    def passed(self) -> bool:
        return not self.checker_mismatches and not self.differential_mismatches


@dataclass
class BenchmarkSuiteResult:
    name: str
    cases: list[BenchmarkCaseResult]

    @property
    def passed(self) -> bool:
        return all(case.passed for case in self.cases)

    def summary(self) -> dict[str, object]:
        return {
            "suite": self.name,
            "passed": self.passed,
            "case_count": len(self.cases),
            "failed_cases": [case.case.name for case in self.cases if not case.passed],
            "total_runs": sum(len(case.runs) for case in self.cases),
        }


__all__ = [
    "RunMetrics",
    "BenchmarkRunSpec",
    "BenchmarkRunResult",
    "DifferentialMismatch",
    "BenchmarkCaseResult",
    "BenchmarkSuiteResult",
]
