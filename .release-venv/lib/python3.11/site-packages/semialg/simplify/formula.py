from __future__ import annotations

import sympy as sp

from .atoms import normalize_atoms
from .boolean import simplify_boolean
from .bounds import simplify_bounds
from .cell_union import cell_union_to_formula
from .equality import simplify_equalities
from .implication import minimize_disj_by_impl

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    sp.PolynomialError,
)


def _canonical_cleanup(expr: sp.Expr) -> sp.Expr:
    """Run inexpensive exact normalizers to a structural fixed point."""

    previous = None
    current = sp.sympify(expr)
    # Equality substitution can expose new comparable bounds, and bound
    # simplification can expose new Boolean identities. Iterate only these
    # inexpensive normalizers to make the result stable without invoking CAD.
    for _ in range(4):
        current = normalize_atoms(current)
        current = simplify_equalities(current)
        current = simplify_bounds(current)
        current = simplify_boolean(current)
        current = normalize_atoms(current)
        if current == previous:
            break
        previous = current
    return current


def simplify_semialgebraic_formula(
    expr: sp.Expr,
    *,
    implication_minimize: bool = True,
    max_dnf_branches: int = 4096,
    max_implication_atoms: int = 8,
    max_implication_vars: int = 3,
    max_implication_degree: int = 4,
) -> sp.Expr:
    """Return a stable canonical simplification of a semialgebraic formula.

    The result is deterministic and idempotent for the supported polynomial
    fragment: polynomial atoms are primitive and consistently oriented,
    repeated zero-set factors are removed, scalar bounds are merged, Boolean
    branches are normalized, and guarded CAD implication checks remove semantic
    redundancies.  This is a stable canonical form, not a claim of globally
    minimum Boolean formula size.
    """

    expr = _canonical_cleanup(expr)
    if implication_minimize:
        try:
            expr = minimize_disj_by_impl(
                expr,
                max_dnf_branches=max_dnf_branches,
                max_atoms=max_implication_atoms,
                max_vars=max_implication_vars,
                max_degree=max_implication_degree,
            )
        except _RECOVERABLE_ERRORS:
            pass
        expr = _canonical_cleanup(expr)
    return expr


def simplify_qe_formula(
    expr: sp.Expr,
    *,
    cell_union=None,
    implication_minimize: bool = True,
    preferred_atoms: tuple[sp.Expr, ...] = (),
    max_dnf_branches: int = 4096,
    max_implication_atoms: int = 8,
    max_implication_vars: int = 3,
    max_implication_degree: int = 4,
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
        except _RECOVERABLE_ERRORS:
            pass
    # Re-introducing original free-variable atoms is exact for conjunctive
    # inputs and gives the simplifier a chance to prefer source polynomials
    # over equivalent root-function cell boundaries.
    if preferred_atoms:
        expr = sp.And(*preferred_atoms, expr)
    return simplify_semialgebraic_formula(
        expr,
        implication_minimize=implication_minimize,
        max_dnf_branches=max_dnf_branches,
        max_implication_atoms=max_implication_atoms,
        max_implication_vars=max_implication_vars,
        max_implication_degree=max_implication_degree,
    )


__all__ = ["simplify_qe_formula", "simplify_semialgebraic_formula"]
