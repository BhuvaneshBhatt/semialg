"""Symbolic validity conditions for canonical geometry objects."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp
from sympy.assumptions.assume import AppliedPredicate

from ._errors import EXACT_OPERATION_ERRORS
from .decision import implies
from .domains import normalize_assumptions


def _predicate_relation(expr: sp.Expr) -> sp.Expr:
    """Translate common SymPy Q predicates into real relational assumptions."""
    if not isinstance(expr, AppliedPredicate) or len(expr.arguments) != 1:
        return expr
    arg = expr.arguments[0]
    predicate = expr.function
    if predicate == sp.Q.positive:
        return arg > 0
    if predicate == sp.Q.nonnegative:
        return arg >= 0
    if predicate == sp.Q.negative:
        return arg < 0
    if predicate == sp.Q.nonpositive:
        return arg <= 0
    if predicate == sp.Q.zero:
        return sp.Eq(arg, 0)
    if predicate == sp.Q.nonzero:
        return sp.Ne(arg, 0)
    return expr


def normalize_geometry_assumptions(
    assumptions: Iterable[sp.Expr] | sp.Expr | None,
) -> tuple[sp.Expr, ...]:
    """Normalize geometry assumptions using the package-wide assumption convention."""
    return tuple(_predicate_relation(item) for item in normalize_assumptions(assumptions))


def _symbols(expressions: Sequence[sp.Expr]) -> tuple[sp.Symbol, ...]:
    return tuple(
        sorted(set().union(*(expr.free_symbols for expr in expressions)), key=sp.default_sort_key)
    )


def _prove(assumptions: tuple[sp.Expr, ...], condition: sp.Expr, *, strategy=None) -> bool:
    condition = sp.simplify(condition)
    if condition is sp.true or condition == sp.true:
        return True
    premise = sp.And(*assumptions) if assumptions else sp.true
    # SymPy's refine/ask path cheaply handles the sign predicates that dominate
    # geometry validity (radii, lengths, principal minors).
    try:
        refined = sp.refine(condition, premise)
        if refined is sp.true or refined == sp.true:
            return True
        if refined is sp.false or refined == sp.false:
            return False
    except (TypeError, ValueError, NotImplementedError):
        pass
    if condition in assumptions:
        return True
    # Only invoke semialg's decision engine when explicitly requested. This
    # keeps ordinary constructor/status queries lightweight and deterministic.
    if strategy is not None:
        variables = _symbols((*assumptions, condition))
        try:
            return bool(implies(premise, condition, variables, strategy=strategy))
        except EXACT_OPERATION_ERRORS:
            return False
    return False


def validity_status(
    conditions: Iterable[sp.Expr],
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    *,
    strategy: str | None = None,
) -> bool | None:
    """Return True, False, or None for validity under the supplied assumptions."""
    required = tuple(sp.simplify(sp.sympify(c)) for c in conditions)
    if not required:
        return True
    asm = normalize_geometry_assumptions(assumptions)
    if all(_prove(asm, condition, strategy=strategy) for condition in required):
        return True
    if any(_prove(asm, sp.Not(condition), strategy=strategy) for condition in required):
        return False
    return None


def validate_geometry_conditions(
    conditions: Iterable[sp.Expr],
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    *,
    label: str = "geometry",
) -> None:
    """Reject a constructor only when its validity conditions are disproved."""
    status = validity_status(conditions, assumptions)
    if status is False:
        raise ValueError(f"{label} validity conditions contradict the supplied assumptions")


__all__ = [
    "normalize_geometry_assumptions",
    "validate_geometry_conditions",
    "validity_status",
]
