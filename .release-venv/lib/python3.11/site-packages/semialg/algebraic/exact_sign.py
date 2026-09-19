"""Exact sign certification for real algebraic expressions."""

from __future__ import annotations

import sympy as sp

from .._zero_testing import certified_zero


def exact_algebraic_sign(expr: object) -> int | None:
    """Return the certified sign of a real algebraic expression when available."""

    value = sp.cancel(sp.sympify(expr))
    if certified_zero(value) is True:
        return 0
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    sign = sp.sign(value)
    if sign in (-1, 0, 1):
        return int(sign)
    try:
        algebraic = sp.polys.numberfields.to_number_field(value)
    except (sp.PolynomialError, TypeError, ValueError, NotImplementedError):
        return None
    if algebraic.is_zero is True:
        return 0
    if algebraic.is_positive is True:
        return 1
    if algebraic.is_negative is True:
        return -1
    sign = sp.sign(algebraic)
    return int(sign) if sign in (-1, 0, 1) else None
