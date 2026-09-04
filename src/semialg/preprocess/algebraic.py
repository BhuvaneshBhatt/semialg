from __future__ import annotations

import sympy as sp
from sympy.core.relational import Relational

from ..structural_keys import ordered_symbols


def _flip_rel(rel: Relational, lhs: sp.Expr) -> sp.Expr:
    if isinstance(rel, sp.Equality):
        return sp.Eq(lhs, 0)
    if isinstance(rel, sp.Unequality):
        return sp.Ne(lhs, 0)
    if isinstance(rel, sp.StrictLessThan):
        return lhs < 0
    if isinstance(rel, sp.LessThan):
        return lhs <= 0
    if isinstance(rel, sp.StrictGreaterThan):
        return lhs > 0
    if isinstance(rel, sp.GreaterThan):
        return lhs >= 0
    return rel


def normalize_polynomial_relation(rel: sp.Expr) -> sp.Expr:
    """Normalize one exact polynomial relation without changing its truth set.

    Equality/disequality atoms are reduced to primitive square-free parts.
    Ordered inequalities remove only certified positive numeric content; a
    negative numeric content is removed while reversing the comparison.
    Parameter-dependent factors are never divided out.
    """

    if not isinstance(rel, Relational):
        return rel
    expr = sp.expand(rel.lhs - rel.rhs)
    symbols = ordered_symbols(expr.free_symbols)
    if not symbols:
        return sp.simplify(rel)
    try:
        poly = sp.Poly(expr, *symbols)
    except (sp.PolynomialError, TypeError, ValueError):
        return rel
    if poly.is_zero:
        return _flip_rel(rel, sp.Integer(0))
    if isinstance(rel, (sp.Equality, sp.Unequality)):
        try:
            _, primitive = poly.primitive()
            normalized = primitive.sqf_part().as_expr()
        except (sp.PolynomialError, TypeError, ValueError):
            normalized = poly.as_expr()
        if sp.Poly(normalized, *symbols).LC().could_extract_minus_sign():
            normalized = -normalized
        return _flip_rel(rel, sp.expand(normalized))
    content, primitive = poly.primitive()
    normalized = sp.expand(primitive.as_expr())
    sign = sp.sign(content)
    reverse = sign == -1
    if sign not in (-1, 1):
        return rel
    if reverse:
        normalized = -normalized
    try:
        leading_negative = bool(sp.Poly(normalized, *symbols).LC().could_extract_minus_sign())
    except (sp.PolynomialError, TypeError, ValueError):
        leading_negative = False
    if leading_negative:
        normalized = -normalized
        reverse = not reverse
    if reverse:
        if isinstance(rel, sp.StrictLessThan):
            return normalized > 0
        if isinstance(rel, sp.LessThan):
            return normalized >= 0
        if isinstance(rel, sp.StrictGreaterThan):
            return normalized < 0
        if isinstance(rel, sp.GreaterThan):
            return normalized <= 0
    return _flip_rel(rel, normalized)


def normalize_polynomial_atoms(expr: sp.Expr) -> sp.Expr:
    """Recursively factor-normalize polynomial atoms in a Boolean formula."""

    if isinstance(expr, Relational):
        return normalize_polynomial_relation(expr)
    if isinstance(expr, sp.And):
        return sp.And(*(normalize_polynomial_atoms(arg) for arg in expr.args))
    if isinstance(expr, sp.Or):
        return sp.Or(*(normalize_polynomial_atoms(arg) for arg in expr.args))
    if isinstance(expr, sp.Not):
        return sp.Not(normalize_polynomial_atoms(expr.args[0]))
    return expr


__all__ = ["normalize_polynomial_atoms", "normalize_polynomial_relation"]
