"""Certified comparison of polynomial models on semialgebraic domains."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..decision import implies
from ..normalization import normalize_formula, normalize_variables
from ..optimization import semialgebraic_maximize, semialgebraic_minimize
from ..optimization_results import OptimizationResult

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


@dataclass(frozen=True)
class PolynomialModelComparisonResult:
    """Exact worst-case discrepancy and dominance information for two models."""

    first_model: sp.Expr
    second_model: sp.Expr
    variables: tuple[sp.Symbol, ...]
    domain: sp.Expr
    difference: sp.Expr
    minimum_difference: OptimizationResult
    maximum_difference: OptimizationResult
    maximum_squared_error: OptimizationResult
    maximum_absolute_error: sp.Expr
    first_le_second: bool
    first_ge_second: bool
    equivalent_on_domain: bool
    counterexamples: Mapping[str, Mapping[sp.Symbol, sp.Expr] | None] = field(default_factory=dict)
    certified: bool = False
    method: str = "exact_polynomial_model_comparison"


def _require_polynomial(expr: sp.Expr, variables: tuple[sp.Symbol, ...], *, label: str) -> None:
    try:
        sp.Poly(expr, *variables)
    except sp.PolynomialError as exc:
        raise ValueError(f"{label} must be polynomial in the predictor variables") from exc


def compare_polynomial_models(
    first: sp.Expr,
    second: sp.Expr,
    variables: Sequence[sp.Symbol | str],
    *,
    domain: FormulaLike | None = None,
    certification: str = "auto",
) -> PolynomialModelComparisonResult:
    """Compare two polynomial models exactly over a semialgebraic domain.

    The comparison reports the exact range endpoints of ``first - second``,
    the maximum absolute discrepancy, dominance in either direction, and exact
    counterexamples to failed dominance claims when available.  Absolute error
    is optimized through the polynomial square ``(first - second)**2`` so the
    core polynomial optimizer remains applicable.
    """

    left = sp.expand(sp.sympify(first))
    right = sp.expand(sp.sympify(second))
    condition = normalize_formula(domain) if domain is not None else sp.true
    vars_ = normalize_variables(
        variables,
        sp.Tuple(left, right, condition),
        append_context_symbols=False,
    )
    if not vars_:
        raise ValueError("at least one predictor variable is required")

    extras = (left.free_symbols | right.free_symbols | condition.free_symbols) - set(vars_)
    if extras:
        names = ", ".join(sorted(symbol.name for symbol in extras))
        raise ValueError(
            "undeclared symbolic parameters are not supported by model comparison: " + names
        )

    _require_polynomial(left, vars_, label="first model")
    _require_polynomial(right, vars_, label="second model")

    difference = sp.expand(left - right)
    min_diff = semialgebraic_minimize(
        difference,
        condition,
        vars_,
        certification=certification,
        return_result=True,
    )
    max_diff = semialgebraic_maximize(
        difference,
        condition,
        vars_,
        certification=certification,
        return_result=True,
    )
    max_sq = semialgebraic_maximize(
        sp.expand(difference**2),
        condition,
        vars_,
        certification=certification,
        return_result=True,
    )
    if not all(isinstance(item, OptimizationResult) for item in (min_diff, max_diff, max_sq)):
        raise TypeError("model comparison requires OptimizationResult values")

    le_check = implies(condition, left <= right, vars_, return_result=True)
    ge_check = implies(condition, left >= right, vars_, return_result=True)
    eq_check = implies(condition, sp.Eq(left, right), vars_, return_result=True)

    max_abs = sp.sqrt(sp.sympify(max_sq.value))
    certified = bool(
        min_diff.certified
        and max_diff.certified
        and max_sq.certified
        and getattr(le_check, "certified", True)
        and getattr(ge_check, "certified", True)
        and getattr(eq_check, "certified", True)
    )
    return PolynomialModelComparisonResult(
        first_model=left,
        second_model=right,
        variables=vars_,
        domain=condition,
        difference=difference,
        minimum_difference=min_diff,
        maximum_difference=max_diff,
        maximum_squared_error=max_sq,
        maximum_absolute_error=sp.simplify(max_abs),
        first_le_second=bool(le_check),
        first_ge_second=bool(ge_check),
        equivalent_on_domain=bool(eq_check),
        counterexamples={
            "first_le_second": getattr(le_check, "counterexample", None)
            if not bool(le_check)
            else None,
            "first_ge_second": getattr(ge_check, "counterexample", None)
            if not bool(ge_check)
            else None,
            "equivalent": getattr(eq_check, "counterexample", None) if not bool(eq_check) else None,
        },
        certified=certified,
    )


__all__ = ["PolynomialModelComparisonResult", "compare_polynomial_models"]
