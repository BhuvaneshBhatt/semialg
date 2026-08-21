"""Certified verification of polynomial Lyapunov functions."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..decision import implies
from ..exact_arithmetic import exact_truth
from ..normalization import normalize_formula, normalize_variables
from ..symbol_resolution import resolve_symbol

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


@dataclass(frozen=True)
class LyapunovVerificationResult:
    """Exact Lyapunov verification result for a polynomial vector field."""

    function: sp.Expr
    variables: tuple[sp.Symbol, ...]
    dynamics: Mapping[sp.Symbol, sp.Expr]
    equilibrium: Mapping[sp.Symbol, sp.Expr]
    domain: sp.Expr
    lie_derivative: sp.Expr
    equilibrium_valid: bool
    equilibrium_in_domain: bool
    positive_definite: bool
    derivative_valid: bool
    valid: bool
    counterexamples: Mapping[str, Mapping[sp.Symbol, sp.Expr] | None] = field(default_factory=dict)
    derivative_strict: bool = True
    method: str = "exact_semialgebraic_lyapunov_verification"
    certified: bool = True


def _normalize_dynamics(
    dynamics: Mapping[sp.Symbol | str, sp.Expr] | Sequence[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    context: sp.Expr,
) -> dict[sp.Symbol, sp.Expr]:
    if isinstance(dynamics, Mapping):
        resolved: dict[sp.Symbol, sp.Expr] = {}
        for raw_var, raw_expr in dynamics.items():
            var = resolve_symbol(raw_var, context=(context,), known_symbols=variables)
            if var not in variables:
                raise ValueError(f"dynamics variable {var!r} is not in the state variable list")
            if var in resolved:
                raise ValueError(f"duplicate dynamics entry for {var!r}")
            resolved[var] = sp.sympify(raw_expr)
    else:
        values = tuple(map(sp.sympify, dynamics))
        if len(values) != len(variables):
            raise ValueError("dynamics sequence must provide one derivative per state variable")
        resolved = dict(zip(variables, values, strict=True))
    missing = [var for var in variables if var not in resolved]
    if missing:
        names = ", ".join(var.name for var in missing)
        raise ValueError(f"dynamics is missing derivatives for: {names}")
    return resolved


def _normalize_equilibrium(
    equilibrium: Mapping[sp.Symbol | str, sp.Expr] | Sequence[sp.Expr] | None,
    variables: tuple[sp.Symbol, ...],
    context: sp.Expr,
) -> dict[sp.Symbol, sp.Expr]:
    if equilibrium is None:
        return {var: sp.Integer(0) for var in variables}
    if isinstance(equilibrium, Mapping):
        resolved: dict[sp.Symbol, sp.Expr] = {}
        for raw_var, raw_value in equilibrium.items():
            var = resolve_symbol(raw_var, context=(context,), known_symbols=variables)
            if var not in variables:
                raise ValueError(f"equilibrium variable {var!r} is not in the state variable list")
            resolved[var] = sp.sympify(raw_value)
        missing = [var for var in variables if var not in resolved]
        if missing:
            names = ", ".join(var.name for var in missing)
            raise ValueError(f"equilibrium is missing values for: {names}")
        return resolved
    values = tuple(map(sp.sympify, equilibrium))
    if len(values) != len(variables):
        raise ValueError("equilibrium sequence must provide one value per state variable")
    return dict(zip(variables, values, strict=True))


def _require_polynomial(expr: sp.Expr, variables: tuple[sp.Symbol, ...], *, label: str) -> None:
    try:
        sp.Poly(expr, *variables)
    except sp.PolynomialError as exc:
        raise ValueError(f"{label} must be polynomial in the state variables") from exc


def verify_lyapunov_function(
    function: sp.Expr,
    dynamics: Mapping[sp.Symbol | str, sp.Expr] | Sequence[sp.Expr],
    variables: Sequence[sp.Symbol | str],
    *,
    domain: FormulaLike | None = None,
    equilibrium: Mapping[sp.Symbol | str, sp.Expr] | Sequence[sp.Expr] | None = None,
    derivative_strict: bool = True,
) -> LyapunovVerificationResult:
    """Certify a proposed polynomial Lyapunov function.

    The function verifies ``V(x*) = 0``, positive definiteness away from the
    equilibrium, and either ``dV/dt < 0`` (default) or ``dV/dt <= 0`` over the
    requested semialgebraic domain.  A strict derivative condition establishes
    the usual asymptotic Lyapunov criterion; the non-strict option verifies the
    weaker Lyapunov-stability derivative condition only.
    """

    candidate = sp.expand(sp.sympify(function))
    dom = normalize_formula(domain) if domain is not None else sp.true
    vars_ = normalize_variables(variables, sp.Tuple(candidate, dom), append_context_symbols=False)
    if not vars_:
        raise ValueError("at least one state variable is required")
    updates = _normalize_dynamics(dynamics, vars_, sp.And(dom, sp.Eq(candidate, candidate)))
    point = _normalize_equilibrium(equilibrium, vars_, sp.And(dom, sp.Eq(candidate, candidate)))

    _require_polynomial(candidate, vars_, label="Lyapunov function")
    for expr in updates.values():
        _require_polynomial(expr, vars_, label="state dynamics")

    expressions = [candidate, dom, *updates.values(), *point.values()]
    extra = set().union(*(expr.free_symbols for expr in expressions)) - set(vars_)
    if extra:
        names = ", ".join(sorted(sym.name for sym in extra))
        raise ValueError(
            "undeclared symbolic parameters are not supported by Lyapunov verification: " + names
        )

    eq_state = all(sp.simplify(updates[var].subs(point)) == 0 for var in vars_)
    eq_value = sp.simplify(candidate.subs(point)) == 0
    try:
        equilibrium_in_domain = exact_truth(dom.subs(point))
    except (TypeError, ValueError, NotImplementedError) as exc:
        raise ValueError(
            "could not verify that the equilibrium lies in the analysis domain"
        ) from exc
    equilibrium_valid = bool(eq_state and eq_value and equilibrium_in_domain)

    distance_sq = sp.expand(sum((var - point[var]) ** 2 for var in vars_))
    away = distance_sq > 0
    positive = implies(sp.And(dom, away), candidate > 0, vars_, return_result=True)

    lie = sp.expand(sum(sp.diff(candidate, var) * updates[var] for var in vars_))
    derivative_goal = lie < 0 if derivative_strict else lie <= 0
    derivative = implies(sp.And(dom, away), derivative_goal, vars_, return_result=True)

    valid = equilibrium_valid and bool(positive) and bool(derivative)
    return LyapunovVerificationResult(
        function=candidate,
        variables=vars_,
        dynamics=updates,
        equilibrium=point,
        domain=dom,
        lie_derivative=lie,
        equilibrium_valid=equilibrium_valid,
        equilibrium_in_domain=equilibrium_in_domain,
        positive_definite=bool(positive),
        derivative_valid=bool(derivative),
        valid=valid,
        counterexamples={
            "positive": getattr(positive, "counterexample", None) if not bool(positive) else None,
            "derivative": getattr(derivative, "counterexample", None)
            if not bool(derivative)
            else None,
        },
        derivative_strict=bool(derivative_strict),
    )


__all__ = ["LyapunovVerificationResult", "verify_lyapunov_function"]
