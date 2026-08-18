from __future__ import annotations

import sympy as sp
from sympy.core.relational import Equality
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import Or as SymOr


def _direct_symbol_solution(eq: Equality) -> tuple[sp.Symbol, sp.Expr] | None:
    lhs, rhs = eq.lhs, eq.rhs
    if isinstance(lhs, sp.Symbol) and lhs not in rhs.free_symbols:
        return lhs, rhs
    if isinstance(rhs, sp.Symbol) and rhs not in lhs.free_symbols:
        return rhs, lhs
    diff = sp.expand(lhs - rhs)
    syms = sorted(diff.free_symbols, key=lambda s: s.name)
    for sym in syms:
        try:
            poly = sp.Poly(diff, sym, domain="EX")
        except Exception:
            continue
        if poly.degree() == 1:
            coeff = poly.coeff_monomial(sym)
            rest = sp.expand(diff - coeff * sym)
            if coeff != 0 and sym not in rest.free_symbols:
                return sym, sp.simplify(-rest / coeff)
    return None


def simplify_equalities(expr: sp.Expr) -> sp.Expr:
    """Substitute simple equalities within each conjunction branch."""

    if isinstance(expr, SymOr):
        return sp.Or(*(simplify_equalities(arg) for arg in expr.args))
    if not isinstance(expr, SymAnd):
        return expr
    args = list(expr.args)
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    equalities: list[sp.Expr] = []
    rest: list[sp.Expr] = []
    for arg in args:
        if isinstance(arg, Equality):
            sol = _direct_symbol_solution(arg)
            if sol is not None:
                sym, value = sol
                substitutions[sym] = value
                equalities.append(sp.Eq(sym, value))
                continue
        rest.append(arg)
    if not substitutions:
        return expr
    simplified_rest = [sp.simplify(arg.subs(substitutions)) for arg in rest]
    return sp.And(*(equalities + simplified_rest))


__all__ = ["simplify_equalities"]
