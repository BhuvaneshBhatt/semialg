"""Very cheap structural reductions used by the public zero oracle."""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from ._cost import within_budget


def literal_number_zero(term: sp.Expr) -> bool | None:
    """Classify finite numeric literals by their represented value.

    Exact SymPy numbers and finite ``Float`` literals are safe to classify
    directly.  A floating literal is inexact as a representation of some
    external mathematical quantity, but its stored value is nevertheless
    exactly zero or nonzero.  Compound expressions containing floats are not
    handled here.
    """
    term = sp.sympify(term)
    if not isinstance(term, sp.Number) or term.has(sp.nan, sp.zoo, sp.oo, -sp.oo):
        return None
    if term.is_finite is not True:
        return None
    zero = term.is_zero
    if zero is not None:
        return bool(zero)
    return bool(term == 0)


def _reduce_node(term: sp.Expr) -> sp.Expr:
    """Apply bounded local arithmetic reductions without general simplification."""
    if not term.args:
        return term
    args = tuple(_reduce_node(arg) for arg in term.args)

    if term.is_Add:
        coeffs = {}
        for arg in args:
            coeff, rest = arg.as_coeff_Mul(rational=True)
            if not coeff.is_Rational:
                return term
            coeffs[rest] = coeffs.get(rest, sp.Integer(0)) + coeff
        parts = []
        for rest, coeff in coeffs.items():
            if coeff == 0:
                continue
            if rest == 1:
                parts.append(coeff)
            elif coeff == 1:
                parts.append(rest)
            else:
                parts.append(sp.Mul(coeff, rest))
        if not parts:
            return sp.Integer(0)
        return sp.Add(*parts)

    if term.is_Mul:
        if any(arg == 0 for arg in args):
            return sp.Integer(0)
        return sp.Mul(*args)

    if term.is_Pow and len(args) == 2:
        base, exp = args
        if exp == 0 and base.is_zero is False:
            return sp.Integer(1)
        if exp == 1:
            return base
        return term

    return term


@lru_cache(maxsize=4096)
def _quick_cached(term: sp.Expr) -> sp.Expr:
    return _reduce_node(term)


def quick_reduce(term: sp.Expr) -> sp.Expr:
    """Apply bounded local arithmetic without generic simplification."""
    term = sp.sympify(term)
    # Preserve a closed compound expression containing approximate operands.
    # Rebuilding it with normal SymPy evaluation can round the operands into a
    # new Float, which is numerical evidence rather than an identity proof.
    if term.has(sp.Float) and not term.free_symbols and not isinstance(term, sp.Number):
        return term
    if not within_budget(term):
        return term
    return _quick_cached(term)
