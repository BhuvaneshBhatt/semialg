"""Certified verification of polynomial barrier certificates."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..decision import implies
from ..normalization import normalize_formula, normalize_variables
from .lyapunov import _normalize_dynamics, _require_polynomial

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


@dataclass(frozen=True)
class BarrierVerificationResult:
    """Exact verification result for a polynomial barrier certificate."""

    barrier: sp.Expr
    variables: tuple[sp.Symbol, ...]
    dynamics: Mapping[sp.Symbol, sp.Expr]
    domain: sp.Expr
    initial_condition: sp.Expr
    unsafe_condition: sp.Expr
    lie_derivative: sp.Expr
    initial_valid: bool
    unsafe_separated: bool
    boundary_valid: bool
    valid: bool
    counterexamples: Mapping[str, Mapping[sp.Symbol, sp.Expr] | None] = field(default_factory=dict)
    derivative_strict: bool = False
    method: str = "exact_semialgebraic_barrier_verification"
    certified: bool = True


def verify_barrier_certificate(
    barrier: sp.Expr,
    dynamics: Mapping[sp.Symbol | str, sp.Expr] | Sequence[sp.Expr],
    variables: Sequence[sp.Symbol | str],
    *,
    initial_condition: FormulaLike,
    unsafe_condition: FormulaLike,
    domain: FormulaLike | None = None,
    derivative_strict: bool = False,
) -> BarrierVerificationResult:
    """Certify a continuous-time polynomial barrier certificate.

    The convention is ``B <= 0`` for the certified safe side.  The verifier
    proves that initial states satisfy ``B <= 0``, unsafe states satisfy
    ``B > 0``, and on the boundary ``B = 0`` the Lie derivative is nonpositive
    (or strictly negative when ``derivative_strict=True``).
    """

    candidate = sp.expand(sp.sympify(barrier))
    dom = normalize_formula(domain) if domain is not None else sp.true
    initial = normalize_formula(initial_condition)
    unsafe = normalize_formula(unsafe_condition)
    context = sp.And(dom, initial, unsafe, sp.Eq(candidate, candidate))
    vars_ = normalize_variables(variables, context, append_context_symbols=False)
    if not vars_:
        raise ValueError("at least one state variable is required")
    updates = _normalize_dynamics(dynamics, vars_, context)

    _require_polynomial(candidate, vars_, label="barrier function")
    for expr in updates.values():
        _require_polynomial(expr, vars_, label="state dynamics")

    expressions = [candidate, dom, initial, unsafe, *updates.values()]
    extra = set().union(*(expr.free_symbols for expr in expressions)) - set(vars_)
    if extra:
        names = ", ".join(sorted(sym.name for sym in extra))
        raise ValueError(
            "undeclared symbolic parameters are not supported by barrier verification: " + names
        )

    initial_check = implies(sp.And(dom, initial), candidate <= 0, vars_, return_result=True)
    unsafe_check = implies(sp.And(dom, unsafe), candidate > 0, vars_, return_result=True)
    lie = sp.expand(sum(sp.diff(candidate, var) * updates[var] for var in vars_))
    derivative_goal = lie < 0 if derivative_strict else lie <= 0
    boundary_check = implies(
        sp.And(dom, sp.Eq(candidate, 0)), derivative_goal, vars_, return_result=True
    )

    valid = bool(initial_check) and bool(unsafe_check) and bool(boundary_check)
    return BarrierVerificationResult(
        barrier=candidate,
        variables=vars_,
        dynamics=updates,
        domain=dom,
        initial_condition=initial,
        unsafe_condition=unsafe,
        lie_derivative=lie,
        initial_valid=bool(initial_check),
        unsafe_separated=bool(unsafe_check),
        boundary_valid=bool(boundary_check),
        valid=valid,
        counterexamples={
            "initial": getattr(initial_check, "counterexample", None)
            if not bool(initial_check)
            else None,
            "unsafe": getattr(unsafe_check, "counterexample", None)
            if not bool(unsafe_check)
            else None,
            "boundary": getattr(boundary_check, "counterexample", None)
            if not bool(boundary_check)
            else None,
        },
        derivative_strict=bool(derivative_strict),
    )


__all__ = ["BarrierVerificationResult", "verify_barrier_certificate"]
