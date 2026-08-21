from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy.core.relational import (
    Equality,
    GreaterThan,
    LessThan,
    Relational,
    StrictGreaterThan,
    StrictLessThan,
)
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import Or as SymOr


@dataclass(frozen=True)
class _Bound:
    value: sp.Expr
    strict: bool
    atom: sp.Expr


def _as_symbol_bound(rel: Relational) -> tuple[sp.Symbol, str, sp.Expr, bool] | None:
    """Return (symbol, side, bound, strict) for simple symbol bounds.

    side is ``lower`` for ``symbol >/>= bound`` and ``upper`` for
    ``symbol </<= bound``. Bounds may contain other symbols but not the bounded
    symbol itself.
    """

    lhs, rhs = rel.lhs, rel.rhs
    strict = isinstance(rel, (StrictLessThan, StrictGreaterThan))
    if isinstance(lhs, sp.Symbol) and lhs not in rhs.free_symbols:
        if isinstance(rel, (StrictLessThan, LessThan)):
            return lhs, "upper", rhs, strict
        if isinstance(rel, (StrictGreaterThan, GreaterThan)):
            return lhs, "lower", rhs, strict
    if isinstance(rhs, sp.Symbol) and rhs not in lhs.free_symbols:
        if isinstance(rel, (StrictLessThan, LessThan)):
            return rhs, "lower", lhs, strict
        if isinstance(rel, (StrictGreaterThan, GreaterThan)):
            return rhs, "upper", lhs, strict
    return None


def _better_lower(left: _Bound, right: _Bound) -> _Bound:
    cmp = sp.simplify(left.value - right.value)
    if cmp.is_positive:
        return left
    if cmp.is_negative:
        return right
    if cmp == 0:
        return left if left.strict or not right.strict else right
    # Unknown symbolic comparison: keep left and let the caller preserve both.
    return left


def _better_upper(left: _Bound, right: _Bound) -> _Bound:
    cmp = sp.simplify(left.value - right.value)
    if cmp.is_negative:
        return left
    if cmp.is_positive:
        return right
    if cmp == 0:
        return left if left.strict or not right.strict else right
    return left


def _bound_atom(sym: sp.Symbol, side: str, bound: _Bound) -> sp.Expr:
    if side == "lower":
        return sym > bound.value if bound.strict else sym >= bound.value
    return sym < bound.value if bound.strict else sym <= bound.value


def _simplify_and_bounds(expr: SymAnd) -> sp.Expr:
    """Simplify a conjunction by combining compatible scalar bounds conservatively."""
    lower: dict[sp.Symbol, _Bound] = {}
    upper: dict[sp.Symbol, _Bound] = {}
    equalities: dict[sp.Symbol, sp.Expr] = {}
    keep: list[sp.Expr] = []
    uncertain_bounds: list[sp.Expr] = []

    for arg in expr.args:
        if isinstance(arg, Equality):
            if isinstance(arg.lhs, sp.Symbol) and arg.lhs not in arg.rhs.free_symbols:
                equalities[arg.lhs] = arg.rhs
                keep.append(arg)
                continue
            if isinstance(arg.rhs, sp.Symbol) and arg.rhs not in arg.lhs.free_symbols:
                equalities[arg.rhs] = arg.lhs
                keep.append(sp.Eq(arg.rhs, arg.lhs))
                continue
        if isinstance(arg, Relational):
            parsed = _as_symbol_bound(arg)
            if parsed is not None:
                sym, side, value, strict = parsed
                candidate = _Bound(value=value, strict=strict, atom=arg)
                if side == "lower":
                    existing = lower.get(sym)
                    if existing is None:
                        lower[sym] = candidate
                    else:
                        best = _better_lower(existing, candidate)
                        if best is existing and sp.simplify(
                            existing.value - candidate.value
                        ) not in (0,):
                            uncertain_bounds.append(arg)
                        lower[sym] = best
                else:
                    existing = upper.get(sym)
                    if existing is None:
                        upper[sym] = candidate
                    else:
                        best = _better_upper(existing, candidate)
                        if best is existing and sp.simplify(
                            existing.value - candidate.value
                        ) not in (0,):
                            uncertain_bounds.append(arg)
                        upper[sym] = best
                continue
        keep.append(arg)

    for sym, value in equalities.items():
        for bound in (lower.get(sym), upper.get(sym)):
            if bound is None:
                continue
            test = sp.simplify(value - bound.value)
            if bound is lower.get(sym):
                if (bound.strict and test.is_positive is False) or (
                    not bound.strict and test.is_nonnegative is False
                ):
                    return sp.false
            else:
                if (bound.strict and test.is_negative is False) or (
                    not bound.strict and test.is_nonpositive is False
                ):
                    return sp.false
        lower.pop(sym, None)
        upper.pop(sym, None)

    atoms = keep + uncertain_bounds
    for sym in sorted(set(lower) | set(upper), key=lambda s: s.name):
        lo = lower.get(sym)
        hi = upper.get(sym)
        if lo is not None and hi is not None:
            cmp = sp.simplify(lo.value - hi.value)
            if cmp.is_positive:
                return sp.false
            if cmp == 0 and (lo.strict or hi.strict):
                return sp.false
        if lo is not None:
            atoms.append(_bound_atom(sym, "lower", lo))
        if hi is not None:
            atoms.append(_bound_atom(sym, "upper", hi))
    return sp.And(*atoms) if atoms else sp.true


def simplify_bounds(expr: sp.Expr) -> sp.Expr:
    """Remove simple redundant symbol bounds in Boolean formulas.

    This intentionally handles a conservative fragment: conjunctions of direct
    symbol lower/upper bounds and simple symbol equalities. More general
    implication-based redundancy should be delegated to CAD later.
    """

    if isinstance(expr, SymOr):
        return sp.Or(*(simplify_bounds(arg) for arg in expr.args))
    if isinstance(expr, SymAnd):
        return _simplify_and_bounds(expr)
    return expr


__all__ = ["simplify_bounds"]
