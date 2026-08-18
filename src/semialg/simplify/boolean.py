from __future__ import annotations

import sympy as sp

from .atoms import normalize_atoms


def _merge_point_with_strict_ray(expr: sp.Expr) -> sp.Expr:
    """Merge simple ``Eq(x, a) | x > a``/``Eq(x, a) | x < a`` patterns.

    ``simplify_logic`` deliberately treats relational atoms propositionally, so
    it does not know that a strict ray plus its endpoint is a closed ray. This
    small semantic cleanup keeps QE output compact without calling scalar
    ``sp.simplify`` on Boolean formulas.
    """

    if not isinstance(expr, sp.Or) or len(expr.args) != 2:
        return expr
    left, right = expr.args
    pairs = ((left, right), (right, left))
    for eq, rel in pairs:
        if not isinstance(eq, sp.Equality):
            continue
        for var, value in ((eq.lhs, eq.rhs), (eq.rhs, eq.lhs)):
            if not isinstance(var, sp.Symbol) or getattr(value, "has", lambda *_: False)(var):
                continue
            if isinstance(rel, sp.StrictGreaterThan):
                if rel.lhs == var and sp.simplify(rel.rhs - value) == 0:
                    return var >= value
                if rel.rhs == var and sp.simplify(rel.lhs - value) == 0:
                    return var <= value
            if isinstance(rel, sp.StrictLessThan):
                if rel.lhs == var and sp.simplify(rel.rhs - value) == 0:
                    return var <= value
                if rel.rhs == var and sp.simplify(rel.lhs - value) == 0:
                    return var >= value
    return expr


def simplify_boolean(expr: sp.Expr) -> sp.Expr:
    """Deterministically simplify a Boolean combination of SymPy relations.

    Keep Boolean formulas in SymPy's logic simplifier. Calling scalar
    ``sp.simplify`` on Boolean ``And``/``Or`` formulas can route through
    arithmetic simplifiers such as ``radsimp`` and trigger deprecation
    warnings about non-expression arguments in ``Mul``.
    """

    expr = normalize_atoms(expr)
    try:
        simplified = sp.simplify_logic(expr, form="dnf")
    except Exception:
        simplified = expr
    return _merge_point_with_strict_ray(simplified)


__all__ = ["simplify_boolean"]
