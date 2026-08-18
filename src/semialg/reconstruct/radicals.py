from __future__ import annotations

import sympy as sp

from .root_functions import root_function_expr


def fiber_root_candidates(poly_expr: sp.Expr, fiber_var: sp.Symbol) -> tuple[sp.Expr, ...]:
    """Return readable linear/quadratic root expressions for a fiber polynomial.

    The order is the ordinary real-root order for the common positive-leading
    quadratic case. For higher degree, or when polynomial conversion fails, the
    caller should fall back to ``root_of``.
    """

    try:
        poly = sp.Poly(sp.expand(poly_expr), fiber_var, domain="EX")
    except Exception:
        return tuple()
    degree = poly.degree()
    if degree <= 0:
        return tuple()
    if degree == 1:
        a, b = poly.all_coeffs()
        if a == 0:
            return tuple()
        return (sp.cancel(-b / a),)
    if degree == 2:
        a, b, c = poly.all_coeffs()
        if a == 0:
            return fiber_root_candidates(sp.expand(b * fiber_var + c), fiber_var)
        disc = sp.factor(b**2 - 4 * a * c)
        sqrt_disc = sp.sqrt(disc)
        low = sp.cancel((-b - sqrt_disc) / (2 * a))
        high = sp.cancel((-b + sqrt_disc) / (2 * a))
        if a.is_negative:
            return (high, low)
        return (low, high)
    return tuple()


def fiber_root_expr(
    poly_expr: sp.Expr | None, fiber_var: sp.Symbol, root_index: int | None
) -> sp.Expr | None:
    """Return a radical/linear root expression, with ``root_of`` fallback."""

    if poly_expr is None or root_index is None:
        return None
    if fiber_var not in sp.sympify(poly_expr).free_symbols:
        return None
    roots = fiber_root_candidates(poly_expr, fiber_var)
    if len(roots) == 1 and root_index >= 0:
        return roots[0]
    if 0 <= root_index < len(roots):
        return roots[root_index]
    # Some lifting paths use one-based root numbers for section labels.
    if 1 <= root_index <= len(roots):
        return roots[root_index - 1]
    return root_function_expr(poly_expr, fiber_var, root_index)


__all__ = ["fiber_root_candidates", "fiber_root_expr"]
