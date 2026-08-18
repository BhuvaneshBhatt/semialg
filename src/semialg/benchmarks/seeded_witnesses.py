from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..instances import find_random_section_wit, sample_free_assignments, sample_modular_points


@dataclass(frozen=True)
class SeededWitnessBatch:
    seed: int
    assignments: tuple[dict[sp.Symbol, object], ...]
    label: str


def gen_seeded_bench_wits(
    variables: Sequence[sp.Symbol],
    *,
    seed: int = 0,
    sample_count: int = 5,
    domain_rules: Mapping[sp.Symbol, Sequence[object]] | None = None,
    modulus: int | None = None,
    label: str = "seeded_free_assignments",
) -> SeededWitnessBatch:
    return SeededWitnessBatch(
        seed,
        tuple(
            sample_free_assignments(
                variables,
                domain_rules=domain_rules,
                sample_count=sample_count,
                modulus=modulus,
                seed=seed,
            )
        ),
        label,
    )


def gen_seeded_points(
    variable_count: int, modulus: int, *, seed: int = 0, sample_count: int = 10
) -> SeededWitnessBatch:
    vars_ = tuple(sp.Symbol(f"x{i}", integer=True) for i in range(variable_count))
    assignments = tuple(
        {var: value for var, value in zip(vars_, pt, strict=True)}
        for pt in sample_modular_points(variable_count, modulus, sample_count, seed=seed)
    )
    return SeededWitnessBatch(seed, assignments, "seeded_modular_points")


def gen_seeded_section_wit(
    formula: sp.Expr, variables: Sequence[sp.Symbol], *, seed: int = 7, attempts: int = 7
) -> SeededWitnessBatch:
    section = find_random_section_wit(formula, variables, seed=seed, attempts=attempts)
    assignments = tuple([section.assignment] if section.assignment is not None else [])
    return SeededWitnessBatch(seed, assignments, "seeded_section_search")


__all__ = [
    "SeededWitnessBatch",
    "gen_seeded_bench_wits",
    "gen_seeded_points",
    "gen_seeded_section_wit",
]
