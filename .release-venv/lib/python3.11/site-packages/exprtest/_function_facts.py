"""Reviewed function facts used by cheap symbolic predicates.

The registry contains a deliberately small set of SymPy-native domain, sign,
and value-class facts.  Rules are admitted only when their semantics are clear
and inexpensive to check.
"""

from __future__ import annotations

import sympy as sp

from ._assumptions import (
    is_integer_value,
    is_negative,
    is_nonnegative,
    is_nonpositive,
    is_nonzero,
    is_positive,
    is_real_value,
)

PROPERTY_FUNCS = frozenset({
    sp.Abs, sp.log, sp.sin, sp.cos, sp.tan, sp.sinc, sp.sinh, sp.cosh, sp.tanh, sp.atan,
    sp.erf, sp.erfc, sp.asin, sp.acos, sp.gamma, sp.factorial, sp.floor,
    sp.ceiling, sp.loggamma, sp.zeta,
})


def _and(*vals: bool | None) -> bool | None:
    if any(v is False for v in vals):
        return False
    return True if vals and all(v is True for v in vals) else None


def _or(*vals: bool | None) -> bool | None:
    if any(v is True for v in vals):
        return True
    return False if vals and all(v is False for v in vals) else None


def _not(value: bool | None) -> bool | None:
    return None if value is None else not value


def _proved(value: bool | None) -> bool | None:
    """Convert a sufficient-condition result into one-sided proof evidence."""
    return True if value is True else None


def _interval(arg, lo, hi, assumptions) -> bool | None:
    return _and(is_nonnegative(arg - lo, assumptions), is_nonnegative(hi - arg, assumptions))


def _domain_unary(term: sp.Expr, assumptions) -> bool | None:
    """Return the reviewed complex-domain condition for reviewed functions."""
    if len(term.args) != 1:
        return None
    x = term.args[0]
    f = term.func
    if f in {sp.Abs, sp.sin, sp.cos, sp.sinh, sp.cosh, sp.erf, sp.erfc, sp.sinc, sp.asin, sp.acos}:
        return True
    if f is sp.log:
        return is_nonzero(x, assumptions)
    if f is sp.atan:
        return _and(is_nonzero(x - sp.I, assumptions), is_nonzero(x + sp.I, assumptions))
    if f is sp.tan:
        return _not(is_integer_value(sp.Rational(1, 2) + x / sp.pi, assumptions))
    if f is sp.tanh:
        return _not(is_integer_value(sp.Rational(1, 2) - sp.I * x / sp.pi, assumptions))
    if f in {sp.gamma, sp.loggamma}:
        noninteger = _not(is_integer_value(x, assumptions))
        return _or(noninteger, is_positive(sp.re(x), assumptions))
    if f is sp.factorial:
        noninteger = _not(is_integer_value(x, assumptions))
        return _or(noninteger, is_nonnegative(sp.re(x), assumptions))
    if f is sp.zeta:
        return is_nonzero(x - 1, assumptions)
    return None


def function_defined(term: sp.Expr, assumptions=True) -> bool | None:
    """Return a reviewed function-domain decision, ignoring child definedness."""
    term = sp.sympify(term)
    if term.func not in PROPERTY_FUNCS:
        return None
    return _domain_unary(term, assumptions)


def function_property(term: sp.Expr, prop: str, assumptions=True) -> bool | None:
    """Evaluate a reviewed function property of a SymPy function call.

    Supported properties are ``integer``, ``real``, ``positive``, ``negative``,
    ``nonnegative``, and ``nonpositive``.  ``None`` means that the function-fact
    registry has no applicable theorem or its condition is undecidable.
    """
    term = sp.sympify(term)
    if term.func not in PROPERTY_FUNCS or len(term.args) != 1:
        return None
    if function_defined(term, assumptions) is False:
        return None
    x = term.args[0]
    f = term.func

    if f is sp.Abs:
        if prop == "integer":
            return _proved(is_integer_value(x, assumptions))
        if prop == "real":
            return True
        if prop == "positive":
            return _proved(is_nonzero(x, assumptions))
        if prop == "nonnegative":
            return True

    if f is sp.log:
        if prop == "real":
            return _proved(is_positive(x, assumptions))
        if prop == "negative":
            return _proved(_and(is_positive(x, assumptions), is_negative(x - 1, assumptions)))
        if prop == "nonpositive":
            return _proved(_and(is_positive(x, assumptions), is_nonpositive(x - 1, assumptions)))
        if prop == "nonnegative":
            return _proved(_and(is_positive(x, assumptions), is_nonnegative(x - 1, assumptions)))
        if prop == "positive":
            return _proved(_and(is_positive(x, assumptions), is_positive(x - 1, assumptions)))

    if f in {sp.sin, sp.cos, sp.sinh} and prop == "real":
        return _proved(is_real_value(x, assumptions))

    if f is sp.cosh:
        if prop == "real":
            return _proved(is_real_value(x, assumptions))
        if prop in {"positive", "nonnegative"}:
            return _proved(is_real_value(x, assumptions))

    if f in {sp.tanh, sp.atan, sp.erf}:
        if prop == "real":
            return _proved(is_real_value(x, assumptions))
        if prop == "negative":
            return _proved(is_negative(x, assumptions))
        if prop == "nonpositive":
            return _proved(is_nonpositive(x, assumptions))
        if prop == "nonnegative":
            return _proved(is_nonnegative(x, assumptions))
        if prop == "positive":
            return _proved(is_positive(x, assumptions))

    if f is sp.erfc:
        if prop == "real":
            return _proved(is_real_value(x, assumptions))
        if prop in {"positive", "nonnegative"}:
            return _proved(is_real_value(x, assumptions))

    if f in {sp.asin, sp.acos}:
        in_dom = _interval(x, -1, 1, assumptions)
        if prop == "real":
            return _proved(in_dom)
        if f is sp.asin:
            if prop == "negative":
                return _proved(_and(in_dom, is_negative(x, assumptions)))
            if prop == "nonpositive":
                return _proved(_and(in_dom, is_nonpositive(x, assumptions)))
            if prop == "nonnegative":
                return _proved(_and(in_dom, is_nonnegative(x, assumptions)))
            if prop == "positive":
                return _proved(_and(in_dom, is_positive(x, assumptions)))
        else:
            if prop == "nonnegative":
                return _proved(in_dom)
            if prop == "positive":
                return _proved(_and(in_dom, is_positive(1 - x, assumptions)))

    if f is sp.gamma:
        if prop == "integer":
            return _proved(_and(is_integer_value(x, assumptions), is_positive(x, assumptions)))
        if prop == "real":
            return _proved(_or(is_positive(x, assumptions), _and(is_real_value(x, assumptions), _not(is_integer_value(x, assumptions)))))
        if prop in {"positive", "nonnegative"}:
            return _proved(is_positive(x, assumptions))

    if f is sp.factorial:
        if prop == "integer":
            return _proved(_and(is_integer_value(x, assumptions), is_nonnegative(x, assumptions)))
        if prop == "real":
            return _proved(_or(is_positive(x + 1, assumptions), _and(is_real_value(x, assumptions), _not(is_integer_value(x, assumptions)))))
        if prop in {"positive", "nonnegative"}:
            return _proved(is_positive(x + 1, assumptions))

    if f in {sp.floor, sp.ceiling}:
        if prop == "integer":
            return _proved(is_real_value(x, assumptions))
        if prop == "real":
            return _proved(is_real_value(x, assumptions))
        if f is sp.floor:
            if prop == "negative": return _proved(is_negative(x, assumptions))
            if prop == "nonpositive": return _proved(is_negative(x - 1, assumptions))
            if prop == "nonnegative": return _proved(is_nonnegative(x, assumptions))
            if prop == "positive": return _proved(is_nonnegative(x - 1, assumptions))
        else:
            if prop == "negative": return _proved(is_nonpositive(x + 1, assumptions))
            if prop == "nonpositive": return _proved(is_nonpositive(x, assumptions))
            if prop == "nonnegative": return _proved(is_positive(x + 1, assumptions))
            if prop == "positive": return _proved(is_positive(x, assumptions))

    if f is sp.loggamma:
        if prop == "real":
            return _proved(is_positive(x, assumptions))
        if prop == "negative":
            return _proved(_and(is_positive(x - 1, assumptions), is_positive(2 - x, assumptions)))
        if prop == "nonpositive":
            return _proved(_and(is_nonnegative(x - 1, assumptions), is_nonnegative(2 - x, assumptions)))
        if prop == "nonnegative":
            return _proved(_or(is_nonnegative(x - 2, assumptions), _and(is_positive(x, assumptions), is_nonnegative(1 - x, assumptions))))
        if prop == "positive":
            return _proved(_or(is_positive(x - 2, assumptions), _and(is_positive(x, assumptions), is_positive(1 - x, assumptions))))

    if f is sp.zeta:
        if prop == "real":
            return _proved(is_real_value(x, assumptions))
        if prop == "negative":
            return _proved(_and(is_positive(1 - x, assumptions), is_positive(x + 2, assumptions)))
        if prop == "nonpositive":
            return _proved(_and(is_positive(1 - x, assumptions), is_nonnegative(x + 2, assumptions)))
        if prop in {"positive", "nonnegative"}:
            return _proved(is_positive(x - 1, assumptions))

    return None


def function_nonzero(term: sp.Expr, assumptions=True) -> bool | None:
    """Return a function-fact nonzero proof from reviewed sign properties."""
    if function_defined(term, assumptions) is not True:
        return None
    pos = function_property(term, "positive", assumptions)
    if pos is True:
        return True
    neg = function_property(term, "negative", assumptions)
    if neg is True:
        return True
    return None


DOMAIN_HAZARD_FUNCS = frozenset({
    sp.log, sp.atan, sp.tan, sp.tanh, sp.gamma, sp.factorial, sp.loggamma, sp.zeta,
})
