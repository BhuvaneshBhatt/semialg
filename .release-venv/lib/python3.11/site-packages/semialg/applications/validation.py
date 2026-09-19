"""Certified validation helpers for symbolic-mathematics results."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..decision import equivalent, implies
from ..normalization import normalize_formula, normalize_variables
from ..optimization import function_range

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


@dataclass(frozen=True)
class ValidationResult:
    """Certified validation outcome with an optional exact counterexample."""

    valid: bool
    claim: str
    counterexample: Mapping[sp.Symbol, sp.Expr] | None = None
    method: str = "exact_semialgebraic_validation"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.valid


def validate_identity(
    lhs: sp.Expr,
    rhs: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] | None = None,
) -> ValidationResult:
    """Certify that two expressions are equal under semialgebraic assumptions."""

    left = sp.sympify(lhs)
    right = sp.sympify(rhs)
    if assumptions is None:
        relation = equivalent(sp.Eq(left, right), sp.true, variables, return_result=True)
        return ValidationResult(
            bool(relation),
            "identity",
            getattr(relation, "counterexample", None),
            getattr(relation, "method", "equivalence"),
        )

    premise = normalize_formula(assumptions)
    vars_ = normalize_variables(variables, sp.Tuple(premise, left, right))
    result = implies(premise, sp.Eq(left, right), vars_, return_result=True)
    return ValidationResult(
        bool(result),
        "identity_under_assumptions",
        getattr(result, "counterexample", None),
        getattr(result, "method", "implication"),
    )


def validate_formula_equivalence(
    original: FormulaLike,
    proposed: FormulaLike,
    variables: Sequence[sp.Symbol | str] | None = None,
) -> ValidationResult:
    """Certify that two formulas define the same real semialgebraic set."""

    result = equivalent(original, proposed, variables, return_result=True)
    return ValidationResult(
        bool(result),
        "formula_equivalence",
        getattr(result, "counterexample", None),
        getattr(result, "method", "equivalence"),
        {"failed_direction": getattr(result, "failed_direction", None)},
    )


def validate_range(
    expression: sp.Expr,
    proposed_range: FormulaLike,
    variables: Sequence[sp.Symbol | str],
    *,
    constraints: FormulaLike | Iterable[FormulaLike] | None = None,
    value_symbol: sp.Symbol | str = "t",
) -> ValidationResult:
    """Certify a proposed exact range formula for an expression."""

    exact = function_range(
        expression,
        constraints,
        variables,
        value_symbol=value_symbol,
        return_result=True,
    )
    result = equivalent(exact.formula, proposed_range, (exact.value_symbol,), return_result=True)
    return ValidationResult(
        bool(result),
        "function_range",
        getattr(result, "counterexample", None),
        getattr(result, "method", "equivalence"),
        {"exact_range": exact.formula},
    )


__all__ = [
    "ValidationResult",
    "validate_formula_equivalence",
    "validate_identity",
    "validate_range",
]
