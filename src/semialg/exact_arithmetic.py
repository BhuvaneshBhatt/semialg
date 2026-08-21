from __future__ import annotations

import sympy as sp
from sympy.core.relational import (
    Equality,
    GreaterThan,
    LessThan,
    StrictGreaterThan,
    StrictLessThan,
    Unequality,
)
from sympy.logic.boolalg import And, BooleanFalse, BooleanTrue, Not, Or


def exact_sign(expr: object) -> int:
    """Return the exact sign of a real algebraic expression.

    The function deliberately refuses undecidable/non-algebraic cases instead
    of turning a fixed-precision numerical approximation into an exact result.
    """

    value = sp.simplify(sp.sympify(expr))
    if value == 0 or value.is_zero is True:
        return 0
    sign = sp.sign(value)
    if sign in (-1, 0, 1):
        return int(sign)
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    try:
        algebraic = sp.polys.numberfields.to_number_field(value).to_root()
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError) as exc:
        raise ValueError(f"could not determine exact sign of {sp.sstr(value)}") from exc
    if algebraic == 0 or algebraic.is_zero is True:
        return 0
    if algebraic.is_positive is True:
        return 1
    if algebraic.is_negative is True:
        return -1
    sign = sp.sign(algebraic)
    if sign in (-1, 0, 1):
        return int(sign)
    raise ValueError(f"could not determine exact sign of {sp.sstr(value)}")


def compare_exact_reals(left: object, right: object) -> int:
    """Compare two exact real algebraic values, returning -1, 0, or 1."""

    lhs = sp.sympify(left)
    rhs = sp.sympify(right)
    if lhs == rhs or sp.simplify(lhs - rhs) == 0:
        return 0
    if lhs is sp.oo or rhs is -sp.oo:
        return 1
    if lhs is -sp.oo or rhs is sp.oo:
        return -1
    return exact_sign(sp.simplify(lhs - rhs))


def exact_truth(expr: object) -> bool:
    """Evaluate a fully specialized Boolean/algebraic relation exactly."""

    value = sp.sympify(expr)
    if value is sp.true or isinstance(value, BooleanTrue):
        return True
    if value is sp.false or isinstance(value, BooleanFalse):
        return False
    if isinstance(value, And):
        return all(exact_truth(arg) for arg in value.args)
    if isinstance(value, Or):
        return any(exact_truth(arg) for arg in value.args)
    if isinstance(value, Not):
        return not exact_truth(value.args[0])
    if isinstance(value, Equality):
        return exact_sign(value.lhs - value.rhs) == 0
    if isinstance(value, Unequality):
        return exact_sign(value.lhs - value.rhs) != 0
    if isinstance(value, StrictLessThan):
        return exact_sign(value.lhs - value.rhs) < 0
    if isinstance(value, LessThan):
        return exact_sign(value.lhs - value.rhs) <= 0
    if isinstance(value, StrictGreaterThan):
        return exact_sign(value.lhs - value.rhs) > 0
    if isinstance(value, GreaterThan):
        return exact_sign(value.lhs - value.rhs) >= 0
    simplified = sp.simplify(value)
    if simplified != value:
        return exact_truth(simplified)
    raise ValueError(f"could not determine exact truth of {sp.sstr(value)}")


__all__ = ["compare_exact_reals", "exact_sign", "exact_truth"]
