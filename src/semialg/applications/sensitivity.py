"""Certified derivative-sign and sensitivity analysis for polynomial models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..decision import implies
from ..normalization import normalize_formula, normalize_variables
from ..optimization import function_range
from ..optimization_results import FunctionRangeResult

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


@dataclass(frozen=True)
class SensitivityDirectionResult:
    """Certified sensitivity information for one predictor variable."""

    variable: sp.Symbol
    derivative: sp.Expr
    range_result: FunctionRangeResult
    nonnegative: bool
    nonpositive: bool
    strictly_positive: bool
    strictly_negative: bool
    constant: bool
    classification: str
    certified: bool = True


@dataclass(frozen=True)
class SensitivityAnalysisResult:
    """Exact coordinate-wise sensitivity summary for a polynomial model."""

    model: sp.Expr
    variables: tuple[sp.Symbol, ...]
    domain: sp.Expr
    directions: Mapping[sp.Symbol, SensitivityDirectionResult]
    method: str = "exact_polynomial_sensitivity"
    certified: bool = True


def _classification(
    *, nonnegative: bool, nonpositive: bool, strictly_positive: bool, strictly_negative: bool
) -> str:
    if nonnegative and nonpositive:
        return "constant"
    if strictly_positive:
        return "strictly_increasing"
    if strictly_negative:
        return "strictly_decreasing"
    if nonnegative:
        return "nondecreasing"
    if nonpositive:
        return "nonincreasing"
    return "mixed"


def analyze_polynomial_sensitivity(
    model: sp.Expr,
    variables: Sequence[sp.Symbol | str],
    *,
    domain: FormulaLike | None = None,
) -> SensitivityAnalysisResult:
    """Certify coordinate-wise derivative signs and derivative ranges.

    The monotonicity labels describe behavior along coordinate-line segments
    contained in the supplied domain.  Arbitrary disconnected domains should
    therefore be interpreted through the derivative-sign fields directly.
    """

    expr = sp.expand(sp.sympify(model))
    condition = normalize_formula(domain) if domain is not None else sp.true
    vars_ = normalize_variables(variables, sp.Tuple(expr, condition), append_context_symbols=False)
    if not vars_:
        raise ValueError("at least one predictor variable is required")
    extra = (expr.free_symbols | condition.free_symbols) - set(vars_)
    if extra:
        names = ", ".join(sorted(sym.name for sym in extra))
        raise ValueError(
            "undeclared symbolic parameters are not supported by sensitivity analysis: " + names
        )
    try:
        sp.Poly(expr, *vars_)
    except sp.PolynomialError as exc:
        raise ValueError("sensitivity model must be polynomial in the predictor variables") from exc

    directions: dict[sp.Symbol, SensitivityDirectionResult] = {}
    for var in vars_:
        derivative = sp.expand(sp.diff(expr, var))
        nonnegative = bool(implies(condition, derivative >= 0, vars_))
        nonpositive = bool(implies(condition, derivative <= 0, vars_))
        strictly_positive = bool(implies(condition, derivative > 0, vars_))
        strictly_negative = bool(implies(condition, derivative < 0, vars_))
        range_result = function_range(derivative, condition, vars_, return_result=True)
        if not isinstance(range_result, FunctionRangeResult):
            raise TypeError("sensitivity analysis requires a FunctionRangeResult")
        constant = bool(nonnegative and nonpositive)
        directions[var] = SensitivityDirectionResult(
            variable=var,
            derivative=derivative,
            range_result=range_result,
            nonnegative=nonnegative,
            nonpositive=nonpositive,
            strictly_positive=strictly_positive,
            strictly_negative=strictly_negative,
            constant=constant,
            classification=_classification(
                nonnegative=nonnegative,
                nonpositive=nonpositive,
                strictly_positive=strictly_positive,
                strictly_negative=strictly_negative,
            ),
        )
    return SensitivityAnalysisResult(expr, vars_, condition, directions)


__all__ = [
    "SensitivityAnalysisResult",
    "SensitivityDirectionResult",
    "analyze_polynomial_sensitivity",
]
