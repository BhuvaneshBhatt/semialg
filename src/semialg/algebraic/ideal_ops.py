"""Small exact ideal operations shared by GTZ implementations."""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from .groebner_utils import compute_groebner_basis


def _expanded_nonzero(generators: Sequence[sp.Expr]) -> tuple[sp.Expr, ...]:
    """Expand each generator once and discard zeros."""

    result = []
    for generator in generators:
        expanded = sp.expand(generator)
        if expanded != 0:
            result.append(expanded)
    return tuple(result)


def canonical_qq_basis_uncached(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol], *, modular: bool | None = None
) -> tuple[sp.Expr, ...]:
    vars_ = tuple(variables)
    exprs = _expanded_nonzero(generators)
    if not exprs:
        return tuple()
    gb = compute_groebner_basis(exprs, vars_, order="grevlex", domain=sp.QQ, modular=modular)
    return tuple(sp.expand(p.as_expr()) for p in gb.polys)


def qq_ideal_equal_uncached(
    left: Sequence[sp.Expr],
    right: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    modular: bool | None = None,
) -> bool:
    return canonical_qq_basis_uncached(
        left, variables, modular=modular
    ) == canonical_qq_basis_uncached(right, variables, modular=modular)


def elimination_ideal_qq(
    generators: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    eliminate: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Return the exact elimination ideal over ``QQ``.

    ``eliminate`` must be a subset of ``variables``. The result is a canonical
    Groebner basis in the retained variable order.
    """
    vars_ = tuple(variables)
    removed = tuple(eliminate)
    if not removed or any(symbol not in vars_ for symbol in removed):
        if not removed:
            return canonical_qq_basis_uncached(generators, vars_)
        raise ValueError("elimination variables must belong to variables")
    retained = tuple(symbol for symbol in vars_ if symbol not in set(removed))
    if not retained:
        raise ValueError("elimination must retain at least one variable")
    exprs = _expanded_nonzero(generators)
    if not exprs:
        return tuple()
    gb = sp.groebner(exprs, *removed, *retained, order="lex", domain=sp.QQ)
    eliminated = tuple(
        sp.expand(poly.as_expr())
        for poly in gb.polys
        if not (set(removed) & poly.as_expr().free_symbols)
    )
    return canonical_qq_basis_uncached(eliminated, retained)


def saturate_ideal_qq(
    generators: Sequence[sp.Expr],
    factor: sp.Expr,
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Return ``<generators> : factor**infinity`` exactly over ``QQ``."""
    vars_ = tuple(variables)
    if not vars_:
        raise ValueError("ideal saturation requires ambient variables")
    factor_expr = sp.expand(sp.sympify(factor))
    if factor_expr == 0:
        raise ValueError("saturation factor must be nonzero")
    if factor_expr.free_symbols - set(vars_):
        raise ValueError("saturation factor contains symbols outside variables")
    t = sp.Dummy("ideal_saturation")
    exprs = _expanded_nonzero(generators)
    gb = sp.groebner((*exprs, 1 - t * factor_expr), t, *vars_, order="lex", domain=sp.QQ)
    saturated = tuple(
        sp.expand(poly.as_expr()) for poly in gb.polys if t not in poly.as_expr().free_symbols
    )
    return canonical_qq_basis_uncached(saturated, vars_)


def ideal_intersection_qq(
    left: Sequence[sp.Expr],
    right: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    canonicalizer=canonical_qq_basis_uncached,
) -> tuple[sp.Expr, ...]:
    vars_ = tuple(variables)
    t = sp.Dummy("ideal_intersection")
    generators = tuple(t * sp.expand(g) for g in left) + tuple(
        (1 - t) * sp.expand(g) for g in right
    )
    gb = sp.groebner(generators, t, *vars_, order="lex", domain=sp.QQ)
    eliminated = tuple(
        sp.expand(p.as_expr()) for p in gb.polys if t not in p.as_expr().free_symbols
    )
    return canonicalizer(eliminated, vars_)


def intersection_all_qq(
    ideals: Sequence[Sequence[sp.Expr]],
    variables: Sequence[sp.Symbol],
    *,
    canonicalizer=canonical_qq_basis_uncached,
) -> tuple[sp.Expr, ...]:
    if not ideals:
        return (sp.Integer(1),)
    current = tuple(ideals[0])
    for ideal in ideals[1:]:
        current = ideal_intersection_qq(current, ideal, variables, canonicalizer=canonicalizer)
    return canonicalizer(current, tuple(variables))


__all__ = [
    "canonical_qq_basis_uncached",
    "elimination_ideal_qq",
    "ideal_intersection_qq",
    "intersection_all_qq",
    "qq_ideal_equal_uncached",
    "saturate_ideal_qq",
]
