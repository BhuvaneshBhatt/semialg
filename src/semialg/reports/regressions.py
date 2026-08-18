from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass

from ..benchmarks import BenchmarkSuiteResult


@dataclass(frozen=True)
class RegressionDelta:
    suite: str
    metric: str
    before: int | float | bool
    after: int | float | bool


def compare_suite_summaries(
    before: BenchmarkSuiteResult, after: BenchmarkSuiteResult
) -> dict[str, object]:
    sb = before.summary()
    sa = after.summary()
    deltas = []
    for key in sorted(set(sb) & set(sa)):
        if sb[key] != sa[key]:
            deltas.append(
                RegressionDelta(suite=after.name, metric=key, before=sb[key], after=sa[key])
            )
    return {
        "suite": after.name,
        "changed": bool(deltas),
        "deltas": [asdict(delta) for delta in deltas],
        "before": sb,
        "after": sa,
    }


def compare_named_runs(
    previous: Mapping[str, BenchmarkSuiteResult], current: Mapping[str, BenchmarkSuiteResult]
) -> dict[str, object]:
    names = sorted(set(previous) & set(current))
    return {name: compare_suite_summaries(previous[name], current[name]) for name in names}
