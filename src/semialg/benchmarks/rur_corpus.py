"""Stable zero-dimensional systems for rational-univariate performance work."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class RURBenchmarkCase:
    name: str
    equations: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]


def rur_benchmark_cases() -> tuple[RURBenchmarkCase, ...]:
    """Return small-to-medium exact systems with varied quotient dimensions."""
    x, y, z = sp.symbols("x y z")
    return (
        RURBenchmarkCase("coupled_quadratic", (x**2 - 2, y - x), (x, y)),
        RURBenchmarkCase("circle_diagonal", (x**2 + y**2 - 1, x - y), (x, y)),
        RURBenchmarkCase(
            "grid_3x3", ((x - 1) * (x - 2) * (x - 3), (y + 1) * (y - 1) * (y - 2)), (x, y)
        ),
        RURBenchmarkCase(
            "grid_4x4",
            ((x - 1) * (x - 2) * (x - 3) * (x - 4), (y + 1) * (y - 1) * (y - 2) * (y - 3)),
            (x, y),
        ),
        RURBenchmarkCase("coupled_three_variable", (x**2 - 2, y**2 - 3, z - x - y), (x, y, z)),
    )


__all__ = ["RURBenchmarkCase", "rur_benchmark_cases"]
