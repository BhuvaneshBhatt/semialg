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
from sympy.polys.polyerrors import CoercionFailed

from .._linear_relations import safe_linear_solution
from ..structural_keys import symbol_identity_key

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


def _canonical_polynomial_residual(expr: sp.Expr, op: str) -> tuple[sp.Expr, str]:
    """Return a primitive, deterministically oriented polynomial residual.

    Semialgebraic atoms often reach the simplifier through different algebraic
    paths (projection, substitution, reconstruction).  Expressions that differ
    only by a nonzero rational scalar should therefore acquire the same
    structural representative.  For inequalities a negative scalar reverses
    the relation; for equality/inequality-to-zero it only changes orientation.
    """

    expr = sp.expand(expr)
    symbols = tuple(sorted(expr.free_symbols, key=symbol_identity_key))
    if not symbols:
        return expr, op
    try:
        poly = sp.Poly(expr, *symbols, domain=sp.QQ)
    except (sp.PolynomialError, TypeError, ValueError, CoercionFailed):
        return expr, op

    # ``clear_denoms`` followed by ``primitive`` produces an integer primitive
    # polynomial independent of rational scaling of the input atom.
    _denom, cleared = poly.clear_denoms(convert=True)
    _content, primitive = cleared.primitive()
    if op in {"=", "!="}:
        # Multiplicity does not affect a zero/nonzero set.  Removing repeated
        # factors makes algebraically identical boundaries canonical even when
        # they were produced with different powers.
        primitive = primitive.sqf_part()
    canon = sp.expand(primitive.as_expr())
    lc = primitive.LC()
    if lc < 0:
        canon = -canon
        op = _reverse_op(op)
    return canon, op


def canonicalize_relation(rel: Relational) -> sp.Expr:
    """Normalize one relational atom to a stable polynomial sign test.

    Polynomial residuals are made primitive over ``QQ`` and oriented by their
    leading coefficient in a deterministic symbol order.  Thus, for example,
    ``2*x - 2*y > 0`` and ``y - x < 0`` canonicalize identically.
    """

    op = _relation_op(rel)
    expr = sp.expand(rel.lhs - rel.rhs)
    expr, op = _canonical_polynomial_residual(expr, op)

    # Preserve the useful human-facing ``x < a`` form for affine univariate
    # atoms after canonical polynomial normalization.
    symbols = tuple(expr.free_symbols)
    if len(symbols) == 1:
        symbol = symbols[0]
        try:
            poly = sp.Poly(expr, symbol)
        except (sp.PolynomialError, TypeError, ValueError):
            poly = None
        if poly is not None and poly.degree() == 1:
            coeff = sp.simplify(poly.coeff_monomial(symbol))
            if coeff.is_number and coeff != 0:
                if coeff.is_negative is True:
                    op = _reverse_op(op)
                if coeff.is_positive is True or coeff.is_negative is True:
                    bound = safe_linear_solution(expr, symbol)
                    if bound is None:
                        return _relation_from_normal(expr, op)
                    if op == "=":
                        return sp.Eq(symbol, bound)
                    if op == "!=":
                        return sp.Ne(symbol, bound)
                    if op == "<":
                        return symbol < bound
                    if op == "<=":
                        return symbol <= bound
                    if op == ">":
                        return symbol > bound
                    if op == ">=":
                        return symbol >= bound

    # Non-polynomial relations still receive the previous inexpensive sign
    # normalization so stable behavior is not restricted to polynomial atoms.
    if not expr.is_polynomial(*symbols):
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


__all__ = ["canonicalize_relation", "normalize_atoms"]
