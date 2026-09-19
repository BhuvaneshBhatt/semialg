"""Cheap exact sign certificates on rectangular domains."""

from __future__ import annotations

from collections.abc import Mapping

import sympy as sp

from ..exact_arithmetic import compare_exact_reals
from ..interval_decomposition import finite_real_roots


def polynomial_sign_on_box(
    expression: object,
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]],
) -> int | None:
    """Certify a constant sign for a separably factored polynomial on a box.

    The return value is ``-1``, ``0``, or ``1`` when the sign is certified and
    ``None`` when this inexpensive proof is insufficient.  Factors depending
    on more than one bounded variable are left to stronger exact
    reasoning.  Boundary zeros are allowed when the interior sign is fixed.
    """
    expr = sp.factor(sp.sympify(expression))
    if expr == 0:
        return 0
    try:
        coefficient, factors = sp.factor_list(expr)
    except (sp.PolynomialError, TypeError, ValueError):
        return None
    try:
        sign = compare_exact_reals(coefficient, 0)
    except (TypeError, ValueError, NotImplementedError):
        return None
    if sign == 0:
        return 0
    sign = 1 if sign > 0 else -1
    bounded = set(bounds)
    for factor, multiplicity in factors:
        variables = factor.free_symbols & bounded
        if not variables:
            try:
                factor_sign = compare_exact_reals(factor, 0)
            except (TypeError, ValueError, NotImplementedError):
                return None
        elif len(variables) == 1:
            variable = next(iter(variables))
            lower, upper = bounds[variable]
            if lower in (-sp.oo, sp.oo) or upper in (-sp.oo, sp.oo):
                return None
            try:
                roots = finite_real_roots(factor, variable)
                interior_root = any(
                    compare_exact_reals(root, lower) > 0 and compare_exact_reals(root, upper) < 0
                    for root in roots
                )
                if interior_root:
                    return None
                midpoint = sp.simplify((lower + upper) / 2)
                factor_sign = compare_exact_reals(factor.subs(variable, midpoint), 0)
            except (
                ArithmeticError,
                TypeError,
                ValueError,
                NotImplementedError,
                sp.PolynomialError,
            ):
                return None
        else:
            return None
        if factor_sign == 0:
            return None
        if multiplicity % 2:
            sign *= 1 if factor_sign > 0 else -1
    return sign


__all__ = ["polynomial_sign_on_box"]
