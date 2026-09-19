"""Exact reference oracles for benchmarking numerical optimizers."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..optimization import semialgebraic_maximize, semialgebraic_minimize
from ..optimization_results import OptimizationResult

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


@dataclass(frozen=True)
class OptimizationBenchmark:
    """Exact reference solution for a numerical optimization benchmark."""

    exact_result: OptimizationResult
    constraints: sp.Expr
    certified: bool

    @property
    def exact_value(self) -> sp.Expr:
        return self.exact_result.value


@dataclass(frozen=True)
class NumericOptimizationCheck:
    """Comparison of a numerical optimizer result with an exact benchmark."""

    benchmark: OptimizationBenchmark
    numeric_value: float
    absolute_error: float
    within_tolerance: bool
    feasible: bool | None = None
    diagnostics: Mapping[str, object] = field(default_factory=dict)


def exact_optimization_benchmark(
    objective: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None,
    variables: Sequence[sp.Symbol | str],
    *,
    kind: str = "min",
    certification: str = "auto",
) -> OptimizationBenchmark:
    """Build a certified exact reference optimum for a numerical benchmark."""

    mode = kind.lower()
    solver = semialgebraic_minimize if mode in {"min", "minimum"} else semialgebraic_maximize
    if mode not in {"min", "minimum", "max", "maximum"}:
        raise ValueError("kind must be 'min' or 'max'")
    result = solver(
        objective,
        constraints,
        variables,
        certification=certification,
        return_result=True,
    )
    if not isinstance(result, OptimizationResult):
        raise TypeError("optimization benchmark requires an OptimizationResult")
    return OptimizationBenchmark(
        result, sp.sympify(constraints if constraints is not None else sp.true), result.certified
    )


def validate_numeric_optimization(
    benchmark: OptimizationBenchmark,
    numeric_value: float,
    *,
    atol: float = 1e-8,
) -> NumericOptimizationCheck:
    """Compare a numerical objective value with a certified exact optimum."""

    if atol < 0:
        raise ValueError("atol must be nonnegative")
    exact_float = float(sp.N(benchmark.exact_value, 30))
    error = abs(float(numeric_value) - exact_float)
    return NumericOptimizationCheck(
        benchmark=benchmark,
        numeric_value=float(numeric_value),
        absolute_error=error,
        within_tolerance=error <= atol,
        diagnostics={"exact_value": benchmark.exact_value},
    )


__all__ = [
    "NumericOptimizationCheck",
    "OptimizationBenchmark",
    "exact_optimization_benchmark",
    "validate_numeric_optimization",
]
