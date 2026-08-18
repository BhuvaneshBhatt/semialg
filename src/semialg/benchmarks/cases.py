from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

from ..formula import ParsedPrenexFormula
from ..parser import parse_quantified_formula


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    formula_text: str
    description: str = ""
    tags: tuple[str, ...] = field(default_factory=tuple)
    sample_points: tuple[Mapping[str, object], ...] = field(default_factory=tuple)
    expected_sample_truth: tuple[bool, ...] = field(default_factory=tuple)

    def parse(self) -> ParsedPrenexFormula:
        return parse_quantified_formula(self.formula_text)


def literature_cases() -> tuple[BenchmarkCase, ...]:
    return (
        BenchmarkCase(
            name="open_unit_disk",
            formula_text="x^2 + y^2 < 1",
            description="Open unit disk CAD benchmark.",
            tags=("literature", "basic", "2d"),
            sample_points=({"x": 0, "y": 0}, {"x": 1, "y": 0}),
            expected_sample_truth=(True, False),
        ),
    )


def nullification_cases() -> tuple[BenchmarkCase, ...]:
    return (
        BenchmarkCase(
            name="specialized_zero_leading_coefficient",
            formula_text="exists y. x*y = 0 and y = 1",
            description="Small nullification-style formula.",
            tags=("nullification",),
        ),
    )


def eq_cons_cases() -> tuple[BenchmarkCase, ...]:
    return (
        BenchmarkCase(
            name="circle_equational_constraint",
            formula_text="exists y. x^2 + y^2 = 1 and y >= 0",
            tags=("ec",),
        ),
    )


def variable_ordering_cases() -> tuple[BenchmarkCase, ...]:
    return (
        BenchmarkCase(
            name="parabola_projection",
            formula_text="exists y. y^2 = x",
            tags=("ordering",),
        ),
    )


def tticad_cases() -> tuple[BenchmarkCase, ...]:
    return (
        BenchmarkCase(
            name="two_formula_family_demo",
            formula_text="(x^2 + y^2 = 1 and x >= 0) or (y = 0 and x > 0)",
            tags=("tticad",),
        ),
    )


__all__ = [
    "BenchmarkCase",
    "literature_cases",
    "nullification_cases",
    "eq_cons_cases",
    "variable_ordering_cases",
    "tticad_cases",
]
