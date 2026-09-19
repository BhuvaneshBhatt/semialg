"""Bounded memoization for immutable exact intermediate objects."""

from __future__ import annotations

from functools import lru_cache

import sympy as sp
from sympy.polys.numberfields import minpoly, to_number_field

from . import _config as cfg

_MP_VAR = sp.Symbol("_tn_m")
_CY_VAR = sp.Symbol("_tn_c")


@lru_cache(maxsize=cfg.EXACT_MINPOLY_CACHE_SIZE)
def _minpoly_expr(term: sp.Expr) -> sp.Expr:
    return sp.Poly(minpoly(term, _MP_VAR), _MP_VAR, domain=sp.QQ).as_expr()


def minpoly_for(term: sp.Expr, var: sp.Symbol) -> sp.Poly:
    """Return a cached exact minimal polynomial written in ``var``."""
    expr = _minpoly_expr(sp.sympify(term)).xreplace({_MP_VAR: var})
    return sp.Poly(expr, var, domain=sp.QQ)


@lru_cache(maxsize=cfg.EXACT_CYCLO_CACHE_SIZE)
def _cyclo_expr(order: int) -> sp.Expr:
    expr = sp.cyclotomic_poly(int(order), _CY_VAR)
    return sp.Poly(expr, _CY_VAR, domain=sp.QQ).as_expr()


def cyclo_for(order: int, var: sp.Symbol) -> sp.Poly:
    """Return a cached cyclotomic polynomial written in ``var``."""
    expr = _cyclo_expr(int(order)).xreplace({_CY_VAR: var})
    return sp.Poly(expr, var, domain=sp.QQ)


@lru_cache(maxsize=cfg.EXACT_ROOT_CACHE_SIZE)
def exact_root(term: sp.Expr):
    """Return a cached exact root representation for an algebraic constant."""
    value = to_number_field(sp.sympify(term))
    if isinstance(value, sp.AlgebraicNumber):
        return value.to_root(radicals=False)
    return value


def clear_exact_cache() -> None:
    """Clear all intermediate exact-object caches."""
    from ._algebraic_bounds import clear_bound_cache
    from ._algebraic_model import clear_model_cache
    from ._cyclotomic import _form_cached
    from ._fast import _quick_cached
    from ._special import _special_cached

    _minpoly_expr.cache_clear()
    _cyclo_expr.cache_clear()
    exact_root.cache_clear()
    _form_cached.cache_clear()
    _quick_cached.cache_clear()
    _special_cached.cache_clear()
    clear_model_cache()
    clear_bound_cache()
