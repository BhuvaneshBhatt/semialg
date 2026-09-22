from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy.core.relational import (
    Equality,
    GreaterThan,
    LessThan,
    Relational,
    StrictGreaterThan,
    StrictLessThan,
)

from ._zero_testing import certified_nonzero, certified_sign

_POLY_ERRORS = (ArithmeticError, TypeError, ValueError, NotImplementedError, sp.PolynomialError)


@dataclass(frozen=True)
class LinearInVariable:
    """Exact degree-one representation ``coefficient*variable + remainder``."""

    variable: sp.Symbol
    coefficient: sp.Expr
    remainder: sp.Expr


@dataclass(frozen=True)
class LinearBound:
    """A globally valid one-sided bound obtained without unsafe division."""

    variable: sp.Symbol
    value: sp.Expr
    side: str
    strict: bool


def linear_in_variable(expr: object, variable: sp.Symbol) -> LinearInVariable | None:
    """Recognize an expression of exact degree one in ``variable``."""

    value = sp.expand(sp.sympify(expr))
    try:
        poly = sp.Poly(value, variable, domain="EX")
    except _POLY_ERRORS:
        return None
    if poly.degree() != 1:
        return None
    coefficient = sp.simplify(poly.coeff_monomial(variable))
    remainder = sp.expand(poly.as_expr() - coefficient * variable)
    if variable in remainder.free_symbols:
        return None
    return LinearInVariable(variable, coefficient, remainder)


def safe_linear_solution(expr: object, variable: sp.Symbol) -> sp.Expr | None:
    """Solve ``expr == 0`` only when the linear coefficient cannot vanish."""

    data = linear_in_variable(expr, variable)
    if data is None or not certified_nonzero(data.coefficient):
        return None
    return sp.cancel(-data.remainder / data.coefficient)


def linear_relation_bound(rel: Relational, variable: sp.Symbol) -> LinearBound | None:
    """Extract a globally valid one-sided bound from a linear inequality.

    The coefficient's sign must be globally certified.  Parameter-dependent
    coefficients with an unresolved sign return ``None`` rather
    than dropping coefficient-zero strata or orienting an inequality wrongly.
    """

    if not isinstance(rel, (StrictLessThan, LessThan, StrictGreaterThan, GreaterThan)):
        return None
    expr = sp.expand(rel.lhs - rel.rhs)
    strict = isinstance(rel, (StrictLessThan, StrictGreaterThan))
    # Normalize to expr <= 0 / expr < 0.
    if isinstance(rel, (StrictGreaterThan, GreaterThan)):
        expr = -expr
    data = linear_in_variable(expr, variable)
    if data is None:
        return None
    sign = certified_sign(data.coefficient)
    if sign not in (-1, 1):
        return None
    value = sp.cancel(-data.remainder / data.coefficient)
    return LinearBound(variable, value, "upper" if sign > 0 else "lower", strict)


def safe_equality_solution(rel: Equality, variable: sp.Symbol) -> sp.Expr | None:
    """Solve one equality for ``variable`` without losing exceptional strata."""

    return safe_linear_solution(sp.expand(rel.lhs - rel.rhs), variable)


__all__ = [
    "LinearBound",
    "LinearInVariable",
    "certified_nonzero",
    "certified_sign",
    "linear_in_variable",
    "linear_relation_bound",
    "safe_equality_solution",
    "safe_linear_solution",
]
