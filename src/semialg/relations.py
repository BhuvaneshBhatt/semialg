"""Low-overhead helpers for normalized polynomial relations."""

from __future__ import annotations

import sympy as sp
from sympy.core.relational import (
    Equality,
    GreaterThan,
    LessThan,
    StrictGreaterThan,
    StrictLessThan,
    Unequality,
)


def split_relation(atom: sp.Expr) -> tuple[sp.Expr, str]:
    """Return ``(lhs - rhs, operator)`` for a SymPy relational atom."""

    if isinstance(atom, Equality):
        return sp.expand(atom.lhs - atom.rhs), "=="
    if isinstance(atom, Unequality):
        return sp.expand(atom.lhs - atom.rhs), "!="
    if isinstance(atom, StrictLessThan):
        return sp.expand(atom.lhs - atom.rhs), "<"
    if isinstance(atom, LessThan):
        return sp.expand(atom.lhs - atom.rhs), "<="
    if isinstance(atom, StrictGreaterThan):
        return sp.expand(atom.lhs - atom.rhs), ">"
    if isinstance(atom, GreaterThan):
        return sp.expand(atom.lhs - atom.rhs), ">="
    raise TypeError(f"expected a relational atom, got {atom!r}")


def make_zero_relation(expr: sp.Expr, operator: str) -> sp.Expr:
    """Construct the relation ``expr operator 0``."""

    if operator == "<":
        return expr < 0
    if operator == "<=":
        return expr <= 0
    if operator == ">":
        return expr > 0
    if operator == ">=":
        return expr >= 0
    if operator == "==":
        return sp.Eq(expr, 0)
    if operator == "!=":
        return sp.Ne(expr, 0)
    raise ValueError(f"unsupported relation operator: {operator!r}")
