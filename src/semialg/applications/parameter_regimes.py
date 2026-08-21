"""Exact parameter-regime analysis for solvability and real-root counts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..conditional import ParameterStratifiedResult
from ..normalization import normalize_formula, normalize_variables
from ..parameters import root_count_conditions, solvability_conditions

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


@dataclass(frozen=True)
class ParameterRegimeResult:
    """Exact piecewise behavior over semialgebraic parameter strata."""

    parameters: tuple[sp.Symbol, ...]
    stratified_result: ParameterStratifiedResult
    quantity: str
    problem: sp.Expr
    variables: tuple[sp.Symbol, ...] = ()
    method: str = "exact_parameter_regime_analysis"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    @property
    def regime_count(self) -> int:
        """Number of nonempty parameter regimes."""

        return self.stratified_result.stratum_count

    @property
    def certified(self) -> bool:
        """Whether every reported regime and the partition are certified."""

        return bool(self.stratified_result.certified)

    def select(self, assignments: Mapping[sp.Symbol | str, object]) -> object:
        """Select the exact regime value at a parameter assignment."""

        return self.stratified_result.select(assignments)


def analyze_parameter_regimes(
    constraints: FormulaLike,
    variables: Sequence[sp.Symbol | str],
    parameters: Sequence[sp.Symbol | str],
) -> ParameterRegimeResult:
    """Partition parameter space by real solvability of a semialgebraic system."""

    expr = normalize_formula(constraints)
    vars_ = normalize_variables(variables, expr, append_context_symbols=False)
    params = normalize_variables(parameters, expr, append_context_symbols=False)
    overlap = set(vars_) & set(params)
    if overlap:
        names = ", ".join(sorted(symbol.name for symbol in overlap))
        raise ValueError(f"variables and parameters must be disjoint: {names}")
    extras = expr.free_symbols - set(vars_) - set(params)
    if extras:
        names = ", ".join(sorted(symbol.name for symbol in extras))
        raise ValueError("all symbolic parameters must be declared: " + names)

    stratified = solvability_conditions(
        expr,
        vars_,
        params,
        return_stratified=True,
    )
    if not isinstance(stratified, ParameterStratifiedResult):
        raise TypeError("parameter regime analysis requires a ParameterStratifiedResult")
    return ParameterRegimeResult(
        parameters=params,
        stratified_result=stratified,
        quantity="solvability",
        problem=expr,
        variables=vars_,
        diagnostics={"regime_values": stratified.values},
    )


def analyze_root_count_regimes(
    polynomial: sp.Expr | sp.Poly,
    variable: sp.Symbol | str,
    parameters: Sequence[sp.Symbol | str] | None = None,
) -> ParameterRegimeResult:
    """Partition parameter space by the number of distinct real polynomial roots."""

    expr = polynomial.as_expr() if isinstance(polynomial, sp.Poly) else sp.sympify(polynomial)
    if isinstance(variable, str):
        matches = tuple(symbol for symbol in expr.free_symbols if symbol.name == variable)
        if len(matches) > 1:
            raise ValueError(f"variable name {variable!r} is ambiguous")
        var = matches[0] if matches else sp.Symbol(variable, real=True)
    else:
        var = variable
    if parameters is None:
        params = tuple(sorted(expr.free_symbols - {var}, key=lambda symbol: symbol.name))
    else:
        params = normalize_variables(parameters, expr, append_context_symbols=False)
    if var in params:
        raise ValueError("root variable and parameters must be disjoint")
    extras = expr.free_symbols - {var} - set(params)
    if extras:
        names = ", ".join(sorted(symbol.name for symbol in extras))
        raise ValueError("all symbolic parameters must be declared: " + names)
    try:
        sp.Poly(expr, var)
    except sp.PolynomialError as exc:
        raise ValueError("root-count regime analysis requires a polynomial") from exc

    stratified = root_count_conditions(
        expr,
        var,
        params,
        return_stratified=True,
    )
    if not isinstance(stratified, ParameterStratifiedResult):
        raise TypeError("root-count regime analysis requires a ParameterStratifiedResult")
    return ParameterRegimeResult(
        parameters=params,
        stratified_result=stratified,
        quantity="real_root_count",
        problem=sp.expand(expr),
        variables=(var,),
        diagnostics={"regime_values": stratified.values},
    )


__all__ = [
    "ParameterRegimeResult",
    "analyze_parameter_regimes",
    "analyze_root_count_regimes",
]
