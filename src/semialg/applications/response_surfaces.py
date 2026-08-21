"""Exact analysis workflows for polynomial response-surface models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..normalization import normalize_formula, normalize_variables
from ..optimization import function_range, semialgebraic_maximize, semialgebraic_minimize
from ..optimization_results import FunctionRangeResult, OptimizationResult

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


@dataclass(frozen=True)
class ResponseSurfaceResult:
    """Exact extrema, range, gradient, and threshold sets for a polynomial model."""

    model: sp.Expr
    variables: tuple[sp.Symbol, ...]
    domain: sp.Expr
    minimum: OptimizationResult
    maximum: OptimizationResult
    range_result: FunctionRangeResult
    gradient: tuple[sp.Expr, ...]
    stationary_condition: sp.Expr
    threshold_regions: Mapping[sp.Expr, sp.Expr] = field(default_factory=dict)
    certified: bool = False
    method: str = "exact_polynomial_response_surface"


def analyze_response_surface(
    model: sp.Expr,
    variables: Sequence[sp.Symbol | str],
    *,
    domain: FormulaLike | None = None,
    thresholds: Sequence[sp.Expr] = (),
    certification: str = "auto",
) -> ResponseSurfaceResult:
    """Analyze a polynomial response surface over a semialgebraic domain.

    The result includes exact global minimum and maximum information, the exact
    image/range formula, the symbolic gradient and stationary condition, and
    exact semialgebraic superlevel sets for requested thresholds.
    """

    expr = sp.expand(sp.sympify(model))
    condition = normalize_formula(domain) if domain is not None else sp.true
    vars_ = normalize_variables(variables, sp.Tuple(expr, condition), append_context_symbols=False)
    if not vars_:
        raise ValueError("at least one predictor variable is required")
    extra_symbols = (expr.free_symbols | condition.free_symbols) - set(vars_)
    if extra_symbols:
        names = ", ".join(sorted(sym.name for sym in extra_symbols))
        raise ValueError(
            "undeclared symbolic parameters are not supported by response-surface analysis: "
            + names
        )
    try:
        sp.Poly(expr, *vars_)
    except sp.PolynomialError as exc:
        raise ValueError(
            "response-surface model must be polynomial in the predictor variables"
        ) from exc

    minimum = semialgebraic_minimize(
        expr, condition, vars_, certification=certification, return_result=True
    )
    maximum = semialgebraic_maximize(
        expr, condition, vars_, certification=certification, return_result=True
    )
    range_result = function_range(expr, condition, vars_, method="auto", return_result=True)
    if not isinstance(minimum, OptimizationResult) or not isinstance(maximum, OptimizationResult):
        raise TypeError("response-surface analysis requires OptimizationResult values")
    if not isinstance(range_result, FunctionRangeResult):
        raise TypeError("response-surface analysis requires a FunctionRangeResult")

    gradient = tuple(sp.diff(expr, var) for var in vars_)
    stationary = sp.And(*(sp.Eq(component, 0) for component in gradient)) if gradient else sp.true
    threshold_regions = {
        sp.sympify(level): sp.And(condition, expr >= sp.sympify(level)) for level in thresholds
    }
    return ResponseSurfaceResult(
        model=expr,
        variables=vars_,
        domain=condition,
        minimum=minimum,
        maximum=maximum,
        range_result=range_result,
        gradient=gradient,
        stationary_condition=stationary,
        threshold_regions=threshold_regions,
        certified=bool(minimum.certified and maximum.certified),
    )


__all__ = ["ResponseSurfaceResult", "analyze_response_surface"]
