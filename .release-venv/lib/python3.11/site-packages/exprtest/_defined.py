"""Cheap definedness, non-pole, and denominator-safety proofs."""

from __future__ import annotations

from collections.abc import MutableMapping
from typing import TYPE_CHECKING

import sympy as sp

from ._assumptions import is_finite, is_nonzero, is_positive, normalize_assumptions
from ._function_facts import function_defined

if TYPE_CHECKING:
    from ._call_memo import OracleMemo


def _memo_key(term: sp.Expr, assumptions) -> tuple:
    return (sp.sympify(term), assumptions)


def quick_defined(term: sp.Expr, assumptions=True) -> bool | None:
    """Return whether ``term`` is defined and finite when cheaply decidable."""
    term = sp.sympify(term)
    assumptions = normalize_assumptions(assumptions)
    if term.has(sp.nan, sp.zoo, sp.oo, -sp.oo):
        return False
    if term.is_Atom:
        return is_finite(term, assumptions)
    if term.is_Add or term.is_Mul:
        vals = tuple(quick_defined(arg, assumptions) for arg in term.args)
        if any(v is False for v in vals):
            return False
        return True if vals and all(v is True for v in vals) else None
    if term.is_Pow:
        base, exp = term.args
        base_def = quick_defined(base, assumptions)
        exp_def = quick_defined(exp, assumptions)
        if base_def is False or exp_def is False:
            return False
        if exp.is_negative is True and is_nonzero(base, assumptions) is False:
            return False
        if base == 0 and exp.is_nonpositive is True:
            return False
        if (
            base_def is True
            and exp_def is True
            and (exp.is_integer is True or is_positive(base, assumptions) is True)
        ):
            return True
        return None
    if term.func is sp.log and len(term.args) == 1:
        arg = term.args[0]
        if quick_defined(arg, assumptions) is False:
            return False
        nz = is_nonzero(arg, assumptions)
        if nz is False:
            return False
        return True if nz is True else None
    if term.func is sp.exp and len(term.args) == 1:
        return quick_defined(term.args[0], assumptions)
    if term.func is sp.gamma and len(term.args) == 1:
        arg = term.args[0]
        if arg.is_integer is True and arg.is_nonpositive is True:
            return False
        if quick_defined(arg, assumptions) is True and (
            arg.is_positive is True or arg.is_integer is False
        ):
            return True
        # Fall through to the function-fact registry for assumption-aware cases.
    domain = function_defined(term, assumptions)
    if domain is not None:
        children = tuple(quick_defined(arg, assumptions) for arg in term.args)
        if any(value is False for value in children):
            return False
        if domain is False:
            return False
        if domain is True and children and all(value is True for value in children):
            return True
    value = is_finite(term, assumptions)
    return value


def _square_nonnegative(term: sp.Expr, assumptions, memo, context) -> tuple[bool, bool]:
    if term.is_Rational:
        return term >= 0, term > 0
    if term.is_Pow and term.exp.is_Integer and term.exp.is_even is True and term.exp.is_positive is True:
        base = term.base
        if base.is_real is True:
            nz = quick_defined_nonzero(base, assumptions, memo, context=context)
            return True, nz is True
    if term.is_Mul:
        coeff, rest = term.as_coeff_Mul()
        if coeff.is_Rational and coeff > 0:
            nn, pos = _square_nonnegative(rest, assumptions, memo, context)
            return nn, pos
    return False, False


def quick_defined_nonzero(term: sp.Expr, assumptions=True, memo: MutableMapping | None = None,
                          *, context: OracleMemo | None = None) -> bool | None:
    """Return True only when finiteness and nonvanishing are cheaply proved."""
    term = sp.sympify(term)
    assumptions = normalize_assumptions(assumptions)
    key = _memo_key(term, assumptions)
    cache = context.defined_cache if context is not None else memo
    if context is not None and cache is None:
        context.defined_cache = {}
        cache = context.defined_cache
    if cache is not None and key in cache:
        return cache[key]
    defined = quick_defined(term, assumptions)
    if defined is False:
        result = False
    else:
        nz = is_nonzero(term, assumptions)
        result = True if defined is True and nz is True else None
        if result is None and term.func is sp.exp and defined is True:
            result = True
        elif result is None and term.is_Mul:
            vals = [quick_defined_nonzero(a, assumptions, memo, context=context) for a in term.args]
            if vals and all(v is True for v in vals):
                result = True
        elif result is None and term.is_Add:
            pieces = [_square_nonnegative(a, assumptions, memo, context) for a in term.args]
            if pieces and all(nn for nn, _ in pieces) and any(pos for _, pos in pieces):
                result = True
        elif result is None and term.is_Pow:
            base, exp = term.args
            if quick_defined_nonzero(base, assumptions, memo, context=context) is True and exp.is_finite is not False:
                result = True
    if cache is not None:
        cache[key] = result
    return result
