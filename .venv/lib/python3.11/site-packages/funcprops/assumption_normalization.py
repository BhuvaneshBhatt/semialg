"""Canonical normalization for public assumptions and domain predicates.

This module is intentionally low-level and side-effect free so every funcprops
query, cache key, and downstream package shares exactly the same normalization.
"""

from __future__ import annotations

import sympy as sp
from sympy import S
from sympy.core.sympify import SympifyError


def _basic(value):
    if value is True or value is None:
        return S.true
    if value is False:
        return S.false
    if isinstance(value, sp.Set):
        return value
    try:
        return sp.sympify(value)
    except (SympifyError, TypeError):
        return value


def _normalize_relation(rel):
    """Cheap deterministic relation normalization for cache keys only."""
    if not isinstance(rel, sp.core.relational.Relational):
        return rel
    lhs, rhs = rel.lhs, rel.rhs
    # Prefer a symbol/expression containing symbols on the left and zero or a
    # simpler constant on the right.  Do not call simplify/QE here.
    if not lhs.free_symbols and rhs.free_symbols:
        return rel.reversed
    return rel


def normalize_assumptions(value):
    """Canonicalize Boolean syntax without expensive logical reasoning."""
    value = _basic(value)
    if value in (S.true, S.false) or isinstance(value, sp.Set):
        return value
    if isinstance(value, (sp.And, sp.Or)):
        head = value.func
        flat = []
        for arg in value.args:
            q = normalize_assumptions(arg)
            if isinstance(q, head):
                flat.extend(q.args)
            else:
                flat.append(q)
        ident = S.true if head is sp.And else S.false
        annih = S.false if head is sp.And else S.true
        if annih in flat:
            return annih
        flat = [q for q in flat if q != ident]
        unique = {q: q for q in flat}
        ordered = tuple(sorted(unique.values(), key=sp.default_sort_key))
        if not ordered:
            return ident
        return head(*ordered, evaluate=False)
    if isinstance(value, sp.Not):
        q = normalize_assumptions(value.args[0])
        return sp.Not(q, evaluate=False)
    return _normalize_relation(value)


def normalize_domain(value):
    if isinstance(value, sp.Set):
        return value
    return normalize_assumptions(value)
