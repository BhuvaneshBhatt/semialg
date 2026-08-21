from __future__ import annotations

import sympy as sp
from sympy.core.relational import (
    Equality,
    GreaterThan,
    LessThan,
    Relational,
    StrictGreaterThan,
    StrictLessThan,
    Unequality,
)
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import BooleanFalse, BooleanTrue
from sympy.logic.boolalg import Not as SymNot
from sympy.logic.boolalg import Or as SymOr

_RELATION_MAP = {
    Equality: "=",
    Unequality: "!=",
    StrictLessThan: "<",
    LessThan: "<=",
    StrictGreaterThan: ">",
    GreaterThan: ">=",
}


def _relation_op(rel: Relational) -> str:
    for cls, op in _RELATION_MAP.items():
        if isinstance(rel, cls):
            return op
    raise TypeError(f"unsupported relation: {rel!r}")


def _reverse_op(op: str) -> str:
    return {"<": ">", "<=": ">=", ">": "<", ">=": "<=", "=": "=", "!=": "!="}[op]


def _relation_from_normal(expr: sp.Expr, op: str) -> sp.Expr:
    expr = sp.expand(expr)
    if expr == 0:
        if op in {"=", "<=", ">="}:
            return sp.true
        return sp.false
    if op == "=":
        return sp.Eq(expr, 0)
    if op == "!=":
        return sp.Ne(expr, 0)
    if op == "<":
        return expr < 0
    if op == "<=":
        return expr <= 0
    if op == ">":
        return expr > 0
    if op == ">=":
        return expr >= 0
    raise ValueError(op)


def canonicalize_relation(rel: Relational) -> sp.Expr:
    """Normalize one relational atom to a deterministic polynomial sign test."""

    op = _relation_op(rel)
    expr = sp.expand(rel.lhs - rel.rhs)
    coeff, primitive = sp.primitive(expr)
    if coeff.is_number and coeff != 0:
        if coeff.could_extract_minus_sign():
            primitive = -primitive
            op = _reverse_op(op)
        expr = primitive
    elif expr.could_extract_minus_sign():
        expr = -expr
        op = _reverse_op(op)
    return _relation_from_normal(expr, op)


def normalize_atoms(expr: sp.Expr) -> sp.Expr:
    """Recursively normalize relational atoms in a SymPy Boolean expression."""

    if expr is True or isinstance(expr, BooleanTrue):
        return sp.true
    if expr is False or isinstance(expr, BooleanFalse):
        return sp.false
    if isinstance(expr, Relational):
        return canonicalize_relation(expr)
    if isinstance(expr, SymAnd):
        return sp.And(*(normalize_atoms(arg) for arg in expr.args))
    if isinstance(expr, SymOr):
        return sp.Or(*(normalize_atoms(arg) for arg in expr.args))
    if isinstance(expr, SymNot):
        return sp.Not(normalize_atoms(expr.args[0]))
    return expr


normalize_relation = canonicalize_relation

__all__ = ["canonicalize_relation", "normalize_atoms", "normalize_relation"]
