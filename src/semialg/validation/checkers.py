from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol

import sympy as sp


@dataclass(frozen=True)
class CheckResult:
    checker_name: str
    available: bool
    formula: sp.Expr | None = None
    truth_value: bool | None = None
    status: str = "unknown"
    diagnostics: tuple[str, ...] = ()


class FormulaChecker(Protocol):
    name: str

    def is_available(self) -> bool: ...

    def eliminate(
        self,
        matrix: sp.Expr,
        variables: Sequence[sp.Symbol],
        quantifiers: Sequence[tuple[str, sp.Symbol]],
    ) -> CheckResult: ...


@dataclass(frozen=True)
class SymPyInequalityChecker:
    """Local checker for unquantified and one-variable inequality fragments."""

    name: str = "sympy-reduce-inequalities"

    def is_available(self) -> bool:
        return True

    def eliminate(
        self,
        matrix: sp.Expr,
        variables: Sequence[sp.Symbol],
        quantifiers: Sequence[tuple[str, sp.Symbol]],
    ) -> CheckResult:
        if quantifiers:
            return CheckResult(self.name, True, status="unsupported")
        try:
            formula = sp.reduce_inequalities(matrix, list(variables))
        except Exception as exc:  # pragma: no cover
            return CheckResult(self.name, True, status="error", diagnostics=(repr(exc),))
        truth = True if formula == sp.true else False if formula == sp.false else None
        return CheckResult(
            self.name, True, formula=formula, truth_value=truth, status="complete-for-fragment"
        )


__all__ = ["CheckResult", "FormulaChecker", "SymPyInequalityChecker"]
