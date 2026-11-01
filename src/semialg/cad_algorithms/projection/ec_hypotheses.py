"""Logical EC propagation and the hypotheses needed between reduced lifts."""

from dataclasses import dataclass
from itertools import combinations

import sympy as sp

from ...algebraic.signs import sign_at_sample
from ..polynomial_utils import normalize_poly


@dataclass(frozen=True)
class ECImplication:
    """An explicit EC or a necessary resultant of two already implied ECs."""

    polynomial: sp.Expr
    parents: tuple[sp.Expr, ...] = ()
    eliminated_variable: sp.Symbol | None = None


@dataclass(frozen=True)
class ECLevelHypotheses:
    """The designated EC, projection operator and required degree conditions."""

    level: int
    equational_constraint: sp.Expr | None
    operator: str
    degree_coefficients: tuple[sp.Expr, ...] = ()
    order_polynomials: tuple[sp.Expr, ...] = ()


def propagated_equations(equations, variables):
    """Derive necessary ECs by a descending, acyclic resultant construction.

    Numeric content and repeated factors may be removed, but individual factors
    of a product are never promoted to necessary equations. Specialization and
    leading-coefficient drops therefore retain all exceptional solution strata.
    """
    variables = tuple(variables)
    if not variables:
        return ()

    def normalized(expression):
        poly = normalize_poly(sp.Poly(expression, *variables))
        return None if poly is None else poly.as_expr()

    known = {}
    for equation in equations:
        expr = normalized(equation)
        if expr is not None:
            known[expr] = ECImplication(expr)
    for level in range(len(variables), 1, -1):
        variable = variables[level - 1]
        active = [
            expr
            for expr in known
            if variable in expr.free_symbols
            and not expr.free_symbols.intersection(variables[level:])
        ]
        for left, right in combinations(active, 2):
            result = normalized(sp.resultant(left, right, variable))
            if result is not None and result not in known:
                known[result] = ECImplication(result, (left, right), variable)
    return tuple(known.values())


def ec_parent_admissible(tower, parent):
    """Whether a parent lies on all designated lower-level EC sections."""
    if parent is None:
        return True
    for plan in tower.metadata.get("ec_levels", ()):
        if plan.level <= parent.level and plan.equational_constraint is not None:
            poly = sp.Poly(plan.equational_constraint, *tower.variables[: parent.level])
            if sign_at_sample(poly, parent.sample) != 0:
                return False
    return True


def replay_ec_projection(tower, formula):
    """Rebuild logical EC implications and the degree/order projection families."""
    from ...formula import equational_constraints, formula_polynomials
    from ..polynomial_utils import polynomial_family_key
    from .reduced import build_reduced_proj_tower

    expected = build_reduced_proj_tower(
        formula_polynomials(formula),
        tower.variables,
        theory=tower.metadata["projection"],
        equational_constraints=equational_constraints(formula),
        recursive_ec=True,
    )
    return (
        tower.metadata.get("ec_levels") == expected.tower.metadata["ec_levels"]
        and tower.metadata.get("ec_implications") == expected.tower.metadata["ec_implications"]
        and len(tower.levels) == len(expected.tower.levels)
        and all(
            left.level == right.level
            and left.variable == right.variable
            and left.entries == right.entries
            and polynomial_family_key(left.polynomials) == polynomial_family_key(right.polynomials)
            for left, right in zip(tower.levels, expected.tower.levels, strict=True)
        )
        and polynomial_family_key(tower.original_polynomials)
        == polynomial_family_key(expected.tower.original_polynomials)
    )
