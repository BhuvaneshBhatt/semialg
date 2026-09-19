"""Certified zero and equality decisions for symbolic validation boundaries."""

from __future__ import annotations

import sympy as sp
from exprtest import zerotest


def certified_zero(expr: object, *, assumptions: object = True) -> bool | None:
    """Return whether ``expr`` is identically zero, or ``None`` if undecidable.

    Both Boolean outcomes come from ``exprtest`` certified mode; callers must
    preserve ``None`` whenever a mathematical decision cannot be certified.
    """

    return zerotest(sp.sympify(expr), assumptions=assumptions, confidence="certified")


def certified_equal(left: object, right: object, *, assumptions: object = True) -> bool | None:
    """Return certified equality of two exact expressions, or ``None``."""

    return certified_zero(sp.sympify(left) - sp.sympify(right), assumptions=assumptions)


def certified_pointwise_zero(expr: object, *, assumptions: object = True) -> bool | None:
    """Return a pointwise zero decision, preserving parameter-dependent uncertainty.

    ``zerotest`` decides algebraic identity.  A nonzero symbolic polynomial such
    as ``a`` is not identically zero, but it can still vanish on the parameter
    stratum ``a = 0``.  This helper therefore promotes a certified non-identity
    result to pointwise nonzero only for symbol-free expressions or when SymPy
    assumptions independently prove nonzeroness.
    """

    value = sp.sympify(expr)
    zero = certified_zero(value, assumptions=assumptions)
    if zero is True:
        return True
    if zero is None:
        return None
    if not value.free_symbols:
        return False
    try:
        nonzero = sp.ask(sp.Q.nonzero(value), assumptions)
    except (TypeError, ValueError, AttributeError):
        nonzero = None
    return False if nonzero is True else None


def certified_constant_sign(expr: object, *, assumptions: object = True) -> int | None:
    """Return the exact sign of a symbol-free real expression when certified.

    This is the shared cheap boundary for constant sign recognition.  Symbolic
    expressions return ``None`` unless their sign is established
    directly by assumptions; callers that need parameter stratification must
    not infer a sign from non-identity.
    """

    value = sp.sympify(expr)
    zero = certified_pointwise_zero(value, assumptions=assumptions)
    if zero is True:
        return 0
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    if value.free_symbols:
        return None
    try:
        sign = sp.sign(sp.simplify(value))
    except (TypeError, ValueError, NotImplementedError):
        return None
    return int(sign) if sign in (-1, 0, 1) else None


def certified_nonzero(expr: object, *, assumptions: object = True) -> bool:
    """Return whether ``expr`` is globally certified to be nonzero.

    Symbol assumptions may prove nonzeroness directly.  Otherwise only a
    symbol-free certified nonzero result is promoted to a global conclusion;
    parameter-dependent non-identity is not enough.
    """

    value = sp.sympify(expr)
    if value.is_zero is False:
        return True
    return certified_pointwise_zero(value, assumptions=assumptions) is False


def require_certified_zero(expr: object, *, assumptions: object = True) -> bool:
    """Return ``True`` only when zero is certified; preserve uncertainty as ``False``.

    This helper is for validation predicates whose caller needs a conservative
    yes/no answer, not for mathematical APIs where ``None`` is meaningful.
    """

    return certified_zero(expr, assumptions=assumptions) is True


__all__ = [
    "certified_constant_sign",
    "certified_equal",
    "certified_nonzero",
    "certified_pointwise_zero",
    "certified_zero",
    "require_certified_zero",
]
