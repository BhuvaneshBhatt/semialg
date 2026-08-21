"""Result models and certification policy for exact optimization."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal

import sympy as sp


@dataclass(frozen=True)
class FunctionRangeResult:
    """Exact range summary for a supported semialgebraic image problem."""

    expression: sp.Expr
    formula: sp.Expr
    value_symbol: sp.Symbol
    variables: tuple[sp.Symbol, ...]
    infimum: sp.Expr | None
    supremum: sp.Expr | None
    minimum_attained: bool | None
    maximum_attained: bool | None
    minimizers: tuple[Mapping[sp.Symbol, sp.Expr], ...] = ()
    maximizers: tuple[Mapping[sp.Symbol, sp.Expr], ...] = ()
    method: str = "qe_image"
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    is_interval: bool | None = None
    interval_count: int | None = None

    @property
    def range_condition(self) -> sp.Expr:
        return self.formula

    @property
    def lower_bound(self) -> sp.Expr | None:
        return self.infimum

    @property
    def upper_bound(self) -> sp.Expr | None:
        return self.supremum

    @property
    def lower_bound_attained(self) -> bool | None:
        return self.minimum_attained

    @property
    def upper_bound_attained(self) -> bool | None:
        return self.maximum_attained


@dataclass(frozen=True)
class OptimizationResult:
    """Exact optimum summary for supported semialgebraic problems."""

    objective: sp.Expr
    variables: tuple[sp.Symbol, ...]
    value: sp.Expr
    points: tuple[Mapping[sp.Symbol, sp.Expr], ...]
    attained: bool
    kind: str
    method: str = "critical_point_enumeration"
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    certified: bool = False

    @property
    def point(self) -> Mapping[sp.Symbol, sp.Expr] | None:
        return self.points[0] if self.points else None


@dataclass(frozen=True)
class ParametricOptimizationResult:
    """Exact first-order characterization of a parameter-dependent optimum."""

    objective: sp.Expr
    constraints: sp.Expr
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    value_symbol: sp.Symbol
    kind: str
    formula: sp.Expr
    quantifiers: tuple[tuple[str, sp.Symbol], ...]
    sample_result: OptimizationResult | None = None
    method: str = "parametric_first_order_optimum_relation"
    certified: bool = True
    quantifier_free: bool = False


@dataclass(frozen=True)
class ParametricFunctionRangeResult:
    """Exact first-order characterization of a parameter-dependent range."""

    expression: sp.Expr
    constraints: sp.Expr
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    value_symbol: sp.Symbol
    formula: sp.Expr
    quantifiers: tuple[tuple[str, sp.Symbol], ...]
    method: str = "parametric_first_order_range_relation"
    certified: bool = True
    quantifier_free: bool = False


@dataclass(frozen=True)
class OptimizationCertificationPolicy:
    """Policy controlling expensive exact global range certification."""

    mode: Literal["auto", "complete", "candidate"] = "auto"
    range_cost_limit: int = 2500
    recursion_limit: int = 4

    def __post_init__(self) -> None:
        if self.mode not in {"auto", "complete", "candidate"}:
            raise ValueError(
                "optimization certification mode must be 'auto', 'complete', or 'candidate'"
            )
        if self.range_cost_limit < 0:
            raise ValueError("range_cost_limit must be nonnegative")
        if self.recursion_limit < 0:
            raise ValueError("recursion_limit must be nonnegative")
