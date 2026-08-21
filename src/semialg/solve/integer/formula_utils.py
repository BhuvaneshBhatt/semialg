"""Formula helpers shared by integer solving strategies."""

from __future__ import annotations

import sympy as sp
from sympy.logic.boolalg import And

_RECOVERABLE_ROOT_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    RuntimeError,
    sp.PolynomialError,
)


def conjuncts(expr: sp.Expr) -> list[sp.Expr]:
    """Return the top-level conjuncts of *expr* without altering their order."""

    return list(expr.args) if isinstance(expr, And) else [expr]


def split_equalities(expr: sp.Expr) -> tuple[list[sp.Expr], list[sp.Expr]]:
    """Split top-level conjuncts into equalities and remaining atoms."""

    atoms = conjuncts(expr)
    equalities = [atom for atom in atoms if isinstance(atom, sp.Equality)]
    others = [atom for atom in atoms if not isinstance(atom, sp.Equality)]
    return equalities, others


def integer_roots_with_completeness(
    poly_expr: sp.Expr, var: sp.Symbol
) -> tuple[list[sp.Expr], bool]:
    """Return integer roots and whether the root list is known complete.

    A failed symbolic root computation is deliberately distinguished from a
    successful proof that no integer roots exist.  Callers that use an empty
    root list to prove infeasibility must inspect the second return value.
    """

    try:
        roots = sp.solveset(sp.Eq(poly_expr, 0), var, domain=sp.S.Integers)
        if isinstance(roots, sp.FiniteSet):
            return list(sorted(roots, key=sp.default_sort_key)), True
        if roots is sp.S.EmptySet or roots == sp.S.EmptySet:
            return [], True
    except _RECOVERABLE_ROOT_ERRORS:
        pass

    try:
        poly = sp.Poly(poly_expr, var)
        roots = poly.all_roots()
    except _RECOVERABLE_ROOT_ERRORS:
        return [], False

    values = {sp.simplify(root) for root in roots}
    integer_values = [value for value in values if value.is_integer is True]
    unknown_values = [value for value in values if value.is_integer is None]
    return list(sorted(integer_values, key=sp.default_sort_key)), not unknown_values


def integer_roots(poly_expr: sp.Expr, var: sp.Symbol) -> list[sp.Expr]:
    """Return exactly identified integer roots of a univariate polynomial.

    This compatibility wrapper omits the completeness flag.  Code that uses an
    empty list as a proof of unsatisfiability should call
    :func:`integer_roots_with_completeness` instead.
    """

    roots, _complete = integer_roots_with_completeness(poly_expr, var)
    return roots


__all__ = [
    "conjuncts",
    "integer_roots",
    "integer_roots_with_completeness",
    "split_equalities",
]
