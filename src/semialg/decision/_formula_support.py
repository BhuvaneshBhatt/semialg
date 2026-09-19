from __future__ import annotations

from collections.abc import Sequence

import sympy as sp
from sympy.logic.boolalg import Boolean

from .._errors import EXACT_OPERATION_ERRORS
from ..formulas.boolean import is_false_expr, is_true_expr
from ..inequality_reduction import reduce_conjunctive_inequalities
from ..simplify.boolean import simplify_boolean
from ..structural_keys import symbol_identity_key
from ._metadata import components_formula, one_dim_components


def safe_simplify_expr(expr: sp.Expr) -> sp.Expr:
    """Simplify an expression without sending Boolean formulas to ``radsimp``."""

    if isinstance(expr, Boolean):
        return simplify_boolean(expr)
    return sp.simplify(expr)


def merge_variables(
    first: Sequence[sp.Symbol],
    second: Sequence[sp.Symbol],
    *formulas: sp.Expr,
) -> tuple[sp.Symbol, ...]:
    """Merge explicit variables and free symbols while preserving order."""

    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for source in (first, second):
        for symbol in source:
            if symbol not in seen:
                out.append(symbol)
                seen.add(symbol)
    for formula in formulas:
        for symbol in sorted(getattr(formula, "free_symbols", set()), key=symbol_identity_key):
            if symbol not in seen:
                out.append(symbol)
                seen.add(symbol)
    return tuple(out)


def fast_parameter_conditions(
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
) -> sp.Expr | None:
    """Return fast parameter conditions for common univariate polynomial atoms."""

    if len(variables) != 1 or not isinstance(expr, sp.core.relational.Relational):
        return None
    variable = variables[0]
    residual = sp.expand(expr.lhs - expr.rhs)
    if not residual.free_symbols <= set(parameters) | {variable}:
        return None
    try:
        polynomial = sp.Poly(residual, variable)
    except EXACT_OPERATION_ERRORS:
        return None
    degree = polynomial.degree()
    if isinstance(expr, sp.Equality):
        if degree == 0:
            return sp.Eq(polynomial.as_expr(), 0)
        if degree == 1:
            linear_coefficient, constant = polynomial.all_coeffs()
            return safe_simplify_expr(
                sp.Or(
                    sp.Ne(linear_coefficient, 0),
                    sp.And(sp.Eq(linear_coefficient, 0), sp.Eq(constant, 0)),
                )
            )
        if degree == 2:
            quadratic_coefficient, linear_coefficient, constant = polynomial.all_coeffs()
            discriminant = sp.expand(linear_coefficient**2 - 4 * quadratic_coefficient * constant)
            quadratic = sp.And(sp.Ne(quadratic_coefficient, 0), discriminant >= 0)
            linear = sp.And(sp.Eq(quadratic_coefficient, 0), sp.Ne(linear_coefficient, 0))
            constant_zero = sp.And(
                sp.Eq(quadratic_coefficient, 0),
                sp.Eq(linear_coefficient, 0),
                sp.Eq(constant, 0),
            )
            return safe_simplify_expr(sp.Or(quadratic, linear, constant_zero))
    return None


def fast_solution_formula(
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, str, bool | None]:
    """Return a conservative exact solution formula for common small solves."""

    if is_true_expr(expr):
        return sp.true, "trivial", True
    if is_false_expr(expr):
        return sp.false, "trivial", False
    if len(variables) == 1:
        components = one_dim_components(expr, variables[0])
        if components is not None:
            reduced = components_formula(components)
            return reduced, "one_dimensional_components", bool(components)
        reduced = reduce_conjunctive_inequalities(expr, variables[0])
        if reduced is not None:
            return reduced, "sympy_reduce_inequalities", not is_false_expr(reduced)
    if len(variables) == 2:
        try:
            from ..implicit_geometry import decompose_cylindrical_formula_to_vertical_bounds_2d

            cells = tuple(decompose_cylindrical_formula_to_vertical_bounds_2d(expr, variables))
            if cells:
                return expr, "vertical_bounds_2d", True
        except EXACT_OPERATION_ERRORS:
            pass
    return expr, "cad", None


__all__ = [
    "fast_parameter_conditions",
    "fast_solution_formula",
    "merge_variables",
    "safe_simplify_expr",
]
