"""Algebraic geometry helpers used by exact optimization."""

from __future__ import annotations

from collections.abc import Sequence
from itertools import combinations

import sympy as sp
from sympy.polys.polyerrors import PolynomialError


def polynomial_locus_dimension(
    equations: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> int | None:
    """Return the exact affine dimension of a polynomial locus when computable."""

    vars_ = tuple(variables)
    if not vars_:
        return 0 if all(sp.expand(eq) == 0 for eq in equations) else -1
    eqs = tuple(sp.expand(eq) for eq in equations if sp.expand(eq) != 0)
    if not eqs:
        return len(vars_)
    try:
        basis = sp.groebner(eqs, *vars_, order="grevlex", domain=sp.QQ)
    except (PolynomialError, ValueError, TypeError):
        return None
    if any(
        poly.as_expr().free_symbols.isdisjoint(set(vars_)) and poly.as_expr() != 0
        for poly in basis.polys
    ):
        return -1
    supports: list[frozenset[int]] = []
    for poly in basis.polys:
        exponents = tuple(poly.LM(order=basis.order).exponents)
        support = frozenset(i for i, exponent in enumerate(exponents) if exponent)
        if not support:
            return -1
        supports.append(support)
    if not supports:
        return len(vars_)
    indices = tuple(range(len(vars_)))
    for size in range(len(vars_) + 1):
        for chosen in combinations(indices, size):
            hit = set(chosen)
            if all(hit.intersection(support) for support in supports):
                return len(vars_) - size
    return 0
