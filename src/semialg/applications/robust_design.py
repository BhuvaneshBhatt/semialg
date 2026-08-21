"""Robust parameter and tolerance analysis for semialgebraic models."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..normalization import normalize_formula, normalize_variables
from ..parameters import solvability_conditions

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


@dataclass(frozen=True)
class RobustParameterResult:
    """Exact parameter conditions for feasible, robust, and violating regimes."""

    constraints: sp.Expr
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    feasible_condition: sp.Expr
    robust_condition: sp.Expr
    violation_condition: sp.Expr
    method: str = "exact_quantifier_elimination"
    diagnostics: Mapping[str, object] = field(default_factory=dict)


def _normalize_parameters(
    parameters: Sequence[sp.Symbol | str],
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Symbol, ...]:
    """Resolve parameter names against symbols already present in the problem."""

    return normalize_variables(
        parameters,
        expr,
        append_context_symbols=False,
    )


def robust_parameter_analysis(
    constraints: FormulaLike | Iterable[FormulaLike],
    operating_variables: Sequence[sp.Symbol | str],
    parameters: Sequence[sp.Symbol | str],
    *,
    operating_domain: FormulaLike | Iterable[FormulaLike] | None = None,
) -> RobustParameterResult:
    """Return exact parameter regions for feasibility and robust satisfaction.

    ``feasible_condition`` characterizes parameter values for which at least
    one operating point satisfies ``constraints``. ``robust_condition``
    characterizes parameter values for which every operating point satisfies
    ``constraints``. ``violation_condition`` characterizes parameter values
    admitting at least one counterexample operating point.
    """

    expr = normalize_formula(constraints)
    domain = normalize_formula(operating_domain) if operating_domain is not None else sp.true
    context = sp.Tuple(expr, domain)
    variables = normalize_variables(
        operating_variables,
        context,
        append_context_symbols=False,
    )
    params = _normalize_parameters(parameters, context, variables)
    overlap = set(variables) & set(params)
    if overlap:
        names = ", ".join(sorted(sym.name for sym in overlap))
        raise ValueError(f"operating variables and parameters must be disjoint: {names}")

    feasible = solvability_conditions(sp.And(domain, expr), variables, params)
    violation = solvability_conditions(sp.And(domain, sp.Not(expr)), variables, params)
    robust = sp.simplify_logic(sp.Not(violation))
    return RobustParameterResult(
        constraints=sp.And(domain, expr),
        variables=variables,
        parameters=params,
        feasible_condition=sp.sympify(feasible),
        robust_condition=sp.sympify(robust),
        violation_condition=sp.sympify(violation),
        diagnostics={
            "quantified_variables": len(variables),
            "parameter_count": len(params),
            "operating_domain": domain,
        },
    )


def robust_parameter_region(
    constraints: FormulaLike | Iterable[FormulaLike],
    operating_variables: Sequence[sp.Symbol | str],
    parameters: Sequence[sp.Symbol | str],
    *,
    quantifier: str = "forall",
    operating_domain: FormulaLike | Iterable[FormulaLike] | None = None,
) -> sp.Expr:
    """Return exact parameter conditions under universal or existential use.

    ``quantifier='forall'`` returns parameter values for which all operating
    points satisfy the constraints. ``quantifier='exists'`` returns parameter
    values for which at least one operating point is feasible.
    """

    result = robust_parameter_analysis(
        constraints,
        operating_variables,
        parameters,
        operating_domain=operating_domain,
    )
    mode = quantifier.lower()
    if mode in {"forall", "all", "robust"}:
        return result.robust_condition
    if mode in {"exists", "any", "feasible"}:
        return result.feasible_condition
    raise ValueError("quantifier must be 'forall' or 'exists'")


__all__ = ["RobustParameterResult", "robust_parameter_analysis", "robust_parameter_region"]
