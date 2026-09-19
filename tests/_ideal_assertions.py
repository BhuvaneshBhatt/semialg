"""Shared exact ideal assertions for algebraic tests."""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp


def canonical_qq_basis(generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol]):
    if not generators:
        return tuple()
    gb = sp.groebner(tuple(generators), *tuple(variables), order="grevlex", domain=sp.QQ)
    return tuple(sp.expand(poly.as_expr()) for poly in gb.polys)


def ideals_equal(left, right, variables) -> bool:
    return canonical_qq_basis(left, variables) == canonical_qq_basis(right, variables)
