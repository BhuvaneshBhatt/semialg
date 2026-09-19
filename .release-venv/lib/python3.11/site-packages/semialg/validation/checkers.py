from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol

import sympy as sp

from ..inequality_reduction import reduce_conjunctive_inequalities


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
        if len(variables) != 1:
            return CheckResult(self.name, True, status="unsupported")
        formula = reduce_conjunctive_inequalities(matrix, variables[0])
        if formula is None:
            return CheckResult(self.name, True, status="unsupported")
        truth = True if formula == sp.true else False if formula == sp.false else None
        return CheckResult(
            self.name, True, formula=formula, truth_value=truth, status="complete-for-fragment"
        )


__all__ = ["CheckResult", "FormulaChecker", "SymPyInequalityChecker"]
