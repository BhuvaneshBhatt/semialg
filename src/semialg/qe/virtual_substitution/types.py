from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import sympy as sp

from ...exceptions import QuantifierEliminationError
from ...status import CoverageStatus, SolverStatus


class VirtualSubstitutionError(QuantifierEliminationError):
    """Raised when a formula is outside the supported quadratic virtual-substitution fragment."""


@dataclass(frozen=True)
class QuadraticVirtualSubstitutionResult:
    """Result of one existential virtual-substitution elimination step."""

    formula: sp.Expr
    eliminated_variable: sp.Symbol
    backend: str = "quadratic-virtual-substitution"
    status: CoverageStatus | str = CoverageStatus.COMPLETE


@dataclass(frozen=True)
class _QuadraticPoint:
    """Root expression ``(constant + radical_sign*sqrt(radical))/denominator``.

    The fields use descriptive names for the algebraic-root parameters so
    substitution code can state the relevant sign conditions directly.
    """

    constant: sp.Expr
    radical_sign: sp.Expr
    radical: sp.Expr
    denominator: sp.Expr


@dataclass(frozen=True)
class VirtualSubstitutionQEResult:
    """Result of a planner-level virtual-substitution QE pass.

    This lightweight result records the quantifier-free formula, the variables
    that were actually eliminated, and whether any quantified variables remain
    for a later backend.
    """

    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    free_variables: tuple[sp.Symbol, ...]
    quantified_variables: tuple[sp.Symbol, ...]
    remaining_quantifiers: tuple[tuple[str, sp.Symbol], ...]
    eliminated_variables: tuple[sp.Symbol, ...]
    backend: str = "quadratic-virtual-substitution-qe"
    status: CoverageStatus | str = CoverageStatus.COMPLETE
    is_sentence: bool = False
    truth_value: bool | None = None
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class VirtualSubstitutionWitnessResult:
    """Witness reconstruction result for a quadratic virtual-substitution pass."""

    instance: Mapping[sp.Symbol, sp.Expr] | None
    eliminated_variables: tuple[sp.Symbol, ...]
    reduced_formula: sp.Expr
    reduced_variables: tuple[sp.Symbol, ...]
    backend: str = "quadratic-virtual-substitution-witness"
    status: SolverStatus | str = "satisfied"
    notes: tuple[str, ...] = ()
