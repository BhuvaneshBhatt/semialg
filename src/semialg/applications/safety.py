"""Certified inductive-invariant checks for discrete polynomial systems."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..decision import implies
from ..normalization import normalize_formula, normalize_variables
from ..symbol_resolution import resolve_symbol

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


@dataclass(frozen=True)
class InvariantVerificationResult:
    """Certified initiation, consecution, and safety checks for an invariant."""

    invariant: sp.Expr
    variables: tuple[sp.Symbol, ...]
    transition: Mapping[sp.Symbol, sp.Expr]
    post_invariant: sp.Expr
    inductive: bool
    initial_valid: bool | None
    safe_valid: bool | None
    valid: bool
    counterexamples: Mapping[str, Mapping[sp.Symbol, sp.Expr] | None] = field(default_factory=dict)
    method: str = "exact_semialgebraic_implication"
    certified: bool = True


def _normalize_transition(
    transition: Mapping[sp.Symbol | str, sp.Expr] | Sequence[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    context: sp.Expr,
) -> dict[sp.Symbol, sp.Expr]:
    if isinstance(transition, Mapping):
        resolved: dict[sp.Symbol, sp.Expr] = {}
        for raw_var, raw_expr in transition.items():
            var = resolve_symbol(raw_var, context=(context,), known_symbols=variables)
            if var not in variables:
                raise ValueError(f"transition variable {var!r} is not in the state variable list")
            if var in resolved:
                raise ValueError(f"duplicate transition for {var!r}")
            resolved[var] = sp.sympify(raw_expr)
    else:
        values = tuple(transition)
        if len(values) != len(variables):
            raise ValueError("transition sequence must provide one update per state variable")
        resolved = dict(zip(variables, map(sp.sympify, values), strict=True))
    missing = [var for var in variables if var not in resolved]
    if missing:
        names = ", ".join(var.name for var in missing)
        raise ValueError(f"transition is missing state updates for: {names}")
    return resolved


def verify_polynomial_invariant(
    invariant: FormulaLike,
    transition: Mapping[sp.Symbol | str, sp.Expr] | Sequence[sp.Expr],
    variables: Sequence[sp.Symbol | str],
    *,
    initial_condition: FormulaLike | None = None,
    unsafe_condition: FormulaLike | None = None,
    domain: FormulaLike | None = None,
) -> InvariantVerificationResult:
    """Certify an inductive safety invariant for a discrete polynomial update.

    The core check proves ``domain & invariant -> invariant(next_state)``.  If
    supplied, ``initial_condition`` must imply the invariant and the invariant
    must exclude ``unsafe_condition``.  Counterexamples are retained whenever
    the decision layer can construct them.
    """

    inv = normalize_formula(invariant)
    dom = normalize_formula(domain) if domain is not None else sp.true
    init = normalize_formula(initial_condition) if initial_condition is not None else None
    unsafe = normalize_formula(unsafe_condition) if unsafe_condition is not None else None
    context = sp.And(
        inv, dom, init if init is not None else sp.true, unsafe if unsafe is not None else sp.true
    )
    vars_ = normalize_variables(variables, context, append_context_symbols=False)
    if not vars_:
        raise ValueError("at least one state variable is required")
    updates = _normalize_transition(transition, vars_, context)

    all_exprs = [inv, dom, *(updates.values())]
    if init is not None:
        all_exprs.append(init)
    if unsafe is not None:
        all_exprs.append(unsafe)
    extra_symbols = set().union(*(expr.free_symbols for expr in all_exprs)) - set(vars_)
    if extra_symbols:
        names = ", ".join(sorted(sym.name for sym in extra_symbols))
        raise ValueError(
            "undeclared symbolic parameters are not supported by invariant verification: " + names
        )

    for expr in updates.values():
        try:
            sp.Poly(expr, *vars_)
        except sp.PolynomialError as exc:
            raise ValueError("state updates must be polynomial in the state variables") from exc

    post = inv.xreplace(updates)
    step = implies(sp.And(dom, inv), post, vars_, return_result=True)
    inductive = bool(step)

    init_result = None
    initial_valid: bool | None = None
    if init is not None:
        init_result = implies(sp.And(dom, init), inv, vars_, return_result=True)
        initial_valid = bool(init_result)

    safe_result = None
    safe_valid: bool | None = None
    if unsafe is not None:
        safe_result = implies(sp.And(dom, inv), sp.Not(unsafe), vars_, return_result=True)
        safe_valid = bool(safe_result)

    valid = inductive and initial_valid is not False and safe_valid is not False
    counterexamples = {
        "inductive": getattr(step, "counterexample", None) if not inductive else None,
        "initial": getattr(init_result, "counterexample", None) if initial_valid is False else None,
        "unsafe": getattr(safe_result, "counterexample", None) if safe_valid is False else None,
    }
    return InvariantVerificationResult(
        invariant=inv,
        variables=vars_,
        transition=updates,
        post_invariant=post,
        inductive=inductive,
        initial_valid=initial_valid,
        safe_valid=safe_valid,
        valid=valid,
        counterexamples=counterexamples,
    )


__all__ = ["InvariantVerificationResult", "verify_polynomial_invariant"]
