from __future__ import annotations

import sympy as sp

from .atoms import normalize_atoms
from .boolean import simplify_boolean
from .bounds import simplify_bounds
from .cell_union import cell_union_to_formula
from .equality import simplify_equalities
from .implication import minimize_disj_by_impl


def simp_semialg_expr(expr: sp.Expr, *, implication_minimize: bool = True) -> sp.Expr:
    expr = normalize_atoms(expr)
    expr = simplify_equalities(expr)
    expr = simplify_bounds(expr)
    expr = simplify_boolean(expr)
    if implication_minimize:
        try:
            expr = minimize_disj_by_impl(expr)
        except Exception:
            pass
        expr = normalize_atoms(expr)
        expr = simplify_equalities(expr)
        expr = simplify_bounds(expr)
        expr = simplify_boolean(expr)
    return expr


def simplify_qe_formula(
    expr: sp.Expr, *, cell_union=None, implication_minimize: bool = True
) -> sp.Expr:
    """Simplify a QE result while preserving CAD-derived semantics.

    When a cell union is available, CAD-cell reconstruction is used as the
    semantic source of truth. A guarded implication minimizer then uses the
    complete CAD backend on small branch pairs to remove multivariate redundant
    atoms and redundant DNF branches.
    """

    if cell_union is not None:
        try:
            expr = cell_union_to_formula(cell_union)
        except Exception:
            pass
    return simp_semialg_expr(expr, implication_minimize=implication_minimize)


__all__ = ["simplify_qe_formula", "simp_semialg_expr"]
