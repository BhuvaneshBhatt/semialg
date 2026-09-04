"""Result records returned by exact reasoning helpers."""

from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp


@dataclass(frozen=True)
class SimplifiedSystem:
    """Result of CAD/QE-backed constraint-system simplification."""

    formula: sp.Expr
    constraints: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]
    inconsistent: bool
    removed_redundant: tuple[sp.Expr, ...] = ()
    substitutions: dict[sp.Symbol, sp.Expr] = field(default_factory=dict)
    method: str = "cad_implication_redundancy"
    diagnostics: dict[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return not self.inconsistent


@dataclass(frozen=True)
class AssumptionSimplificationResult:
    """Structured result for assumption-aware expression simplification."""

    expression: sp.Expr
    original: sp.Expr
    assumptions: sp.Expr
    variables: tuple[sp.Symbol, ...]
    conditions: tuple[sp.Expr, ...] = ()
    rewrites: tuple[str, ...] = ()
    diagnostics: dict[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.expression is not sp.nan


@dataclass(frozen=True)
class SignProofResult:
    """Structured result for proving the sign of a real expression."""

    proven: bool
    relation: str
    expression: sp.Expr
    assumptions: sp.Expr
    variables: tuple[sp.Symbol, ...]
    formula: sp.Expr
    counterexample: dict[sp.Symbol, sp.Expr] | None = None
    method: str = "unknown"
    certificate: object | None = None
    diagnostics: dict[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.proven


@dataclass(frozen=True)
class SignClassificationResult:
    """Canonical exact sign classification of an expression on a region."""

    classification: str
    expression: sp.Expr
    assumptions: sp.Expr
    variables: tuple[sp.Symbol, ...]
    proofs: dict[str, SignProofResult] = field(default_factory=dict)
    diagnostics: dict[str, object] = field(default_factory=dict)

    @property
    def certified(self) -> bool:
        return self.classification != "unknown"


__all__ = [
    "SimplifiedSystem",
    "AssumptionSimplificationResult",
    "SignProofResult",
    "SignClassificationResult",
]
