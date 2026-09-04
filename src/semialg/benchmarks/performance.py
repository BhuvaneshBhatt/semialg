"""Deterministic performance-regression probes for semialg core algorithms."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from time import perf_counter

import sympy as sp


@dataclass(frozen=True)
class PerformanceProbe:
    name: str
    seconds: float
    metrics: dict[str, int | float | str]


def run_core_performance_probes() -> tuple[PerformanceProbe, ...]:
    """Run small stable probes used by the dedicated performance workflow.

    Wall times are reported but deterministic work counters are the primary
    regression contract.  Ordinary unit tests should assert the latter rather
    than brittle machine-dependent timing thresholds.
    """

    from ..optimization import semialgebraic_minimize
    from ..simplify.implication import (
        clear_implication_minimization_stats,
        implication_minimization_stats,
        minimize_disj_by_impl,
    )

    x, y = sp.symbols("x y", real=True)
    probes: list[PerformanceProbe] = []

    formula = sp.Or(
        sp.And(x > 0, y > 0),
        sp.And(x > 0, y <= 0),
        evaluate=False,
    )
    clear_implication_minimization_stats()
    start = perf_counter()
    simplified = minimize_disj_by_impl(formula)
    elapsed = perf_counter() - start
    stats = implication_minimization_stats()
    probes.append(
        PerformanceProbe(
            "shared_cad_boolean_reconstruction",
            elapsed,
            {
                "shared_cad_builds": stats.shared_cad_builds,
                "pairwise_fallbacks": stats.pairwise_fallbacks,
                "cover_selected": stats.cover_selected,
                "result": sp.sstr(simplified),
            },
        )
    )

    start = perf_counter()
    result = semialgebraic_minimize(
        x + 2 * y,
        sp.And(x >= 0, x <= 2, y >= -1, y <= 3),
        (x, y),
        return_result=True,
    )
    elapsed = perf_counter() - start
    probes.append(
        PerformanceProbe(
            "linear_box_optimization",
            elapsed,
            {
                "method": result.method,
                "candidate_count": int(result.diagnostics.get("candidate_count", 0)),
            },
        )
    )

    start = perf_counter()
    separable = semialgebraic_minimize(
        (x - 1) ** 2 + (y + 2) ** 2,
        sp.And(x >= 0, x <= 4, y >= -3, y <= 1),
        (x, y),
        return_result=True,
    )
    elapsed = perf_counter() - start
    probes.append(
        PerformanceProbe(
            "separable_box_optimization",
            elapsed,
            {"method": separable.method, "component_count": len(separable.variables)},
        )
    )
    return tuple(probes)


def probes_as_dict() -> list[dict[str, object]]:
    return [asdict(probe) for probe in run_core_performance_probes()]


__all__ = ["PerformanceProbe", "run_core_performance_probes", "probes_as_dict"]
