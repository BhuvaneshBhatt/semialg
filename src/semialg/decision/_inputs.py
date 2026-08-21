from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp
from sympy.logic.boolalg import Boolean
from sympy.polys.polyerrors import PolynomialError

from ..normalization import normalize_formula, normalize_problem_variables

FormulaLike = sp.Expr | Boolean | bool


def as_real_symbol(var: sp.Symbol | str) -> sp.Symbol:
    return sp.Symbol(var, real=True) if isinstance(var, str) else var


def normalize_decision_variables(
    variables: Sequence[sp.Symbol | str] | None,
    formula: sp.Expr,
) -> tuple[sp.Symbol, ...]:
    return normalize_problem_variables(variables, formula)


def normalize_symbols(symbols: Sequence[sp.Symbol | str] | None) -> tuple[sp.Symbol, ...]:
    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for item in symbols or ():
        sym = as_real_symbol(item)
        if sym not in seen:
            out.append(sym)
            seen.add(sym)
    return tuple(out)


def normalize_solve_variables(
    variables: Sequence[sp.Symbol | str] | None,
    formula: sp.Expr,
    parameters: Sequence[sp.Symbol] = (),
) -> tuple[sp.Symbol, ...]:
    params = set(parameters)
    if variables is not None:
        return tuple(sym for sym in normalize_symbols(variables) if sym not in params)
    return tuple(sorted(formula.free_symbols - params, key=lambda item: item.name))


def prepare_solve_inputs(
    constraints: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None,
    parameters: Sequence[sp.Symbol | str] | None,
    *,
    domain: str,
    method: str,
    variable_order: Sequence[sp.Symbol | str] | None,
    projection_order: Sequence[sp.Symbol | str] | None,
    normalize_domains: bool,
) -> tuple[sp.Expr, sp.Expr, tuple[sp.Symbol, ...], tuple[sp.Symbol, ...], str, object | None]:
    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError("solve_semialgebraic currently supports only the real domain")
    original = normalize_formula(constraints)
    params = normalize_symbols(parameters)
    solve_vars = normalize_solve_variables(variables, original, params)
    param_set = set(params)
    if variable_order is not None:
        ordered = tuple(sym for sym in normalize_symbols(variable_order) if sym not in param_set)
        ordered_set = set(ordered)
        solve_vars = ordered + tuple(sym for sym in solve_vars if sym not in ordered_set)
    method_key = method.lower().replace("-", "_")
    allowed = {"auto", "interval", "linear", "rur", "cad", "qe", "cylindrical", "sampling"}
    if method_key not in allowed:
        raise ValueError(f"unsupported solve_semialgebraic method: {method!r}")
    if projection_order is not None and variable_order is None:
        ordered = tuple(sym for sym in normalize_symbols(projection_order) if sym not in param_set)
        ordered_set = set(ordered)
        solve_vars = ordered + tuple(sym for sym in solve_vars if sym not in ordered_set)

    normalized = original
    domain_info = None
    if normalize_domains:
        try:
            from ..domain_solve import normalize_domain_sensitive_constraints

            domain_info = normalize_domain_sensitive_constraints(original, solve_vars)
            normalized = domain_info.formula
        except (TypeError, ValueError, ArithmeticError, NotImplementedError, PolynomialError):
            normalized = original
    solve_vars = normalize_solve_variables(tuple(solve_vars), normalized, params)
    return original, normalized, params, solve_vars, method_key, domain_info


__all__ = [
    "as_real_symbol",
    "normalize_decision_variables",
    "normalize_solve_variables",
    "normalize_symbols",
    "prepare_solve_inputs",
]
