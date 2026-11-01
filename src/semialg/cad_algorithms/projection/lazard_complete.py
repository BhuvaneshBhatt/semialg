"""The full Lazard projection operator, independent of EC reduction."""

from itertools import combinations

import sympy as sp

from ..performance_cache import PROJECTION_TOWERS
from ..polynomial_utils import (
    _cad_input_polynomials,
    discriminant_expression,
    lower_polynomial,
    polynomial_family_key,
    polynomial_key,
    resultant_expression,
    univariate_leading_coefficient,
)
from .collins import ProjectionLevel, ProjectionPolynomial, ProjectionTower, squarefree_basis


def lazard_projection_entries(polys, variable, lower_variables):
    """Project contents, leading/trailing coefficients, discriminants and resultants."""
    lower_variables = tuple(lower_variables)
    active = []
    claims = []

    def add(expression, operation, parents):
        if not lower_variables:
            return
        polynomial = lower_polynomial(expression, lower_variables)
        if polynomial is not None:
            for factor in squarefree_basis((polynomial,)):
                claims.append(
                    ProjectionPolynomial(
                        factor,
                        len(lower_variables),
                        operation,
                        tuple(parents),
                        variable,
                        sp.expand(expression),
                    )
                )

    for poly in squarefree_basis(polys):
        key = polynomial_key(poly)
        if poly.degree(variable) <= 0:
            add(poly.as_expr(), "inactive", (key,))
            continue
        active.append(poly)
        add(univariate_leading_coefficient(poly, variable), "leading-coefficient", (key,))
        add(poly.as_expr().subs(variable, 0), "trailing-coefficient", (key,))
        if poly.degree(variable) > 1:
            add(discriminant_expression(poly, variable), "discriminant", (key,))
    for left, right in combinations(active, 2):
        add(
            resultant_expression(left, right, variable),
            "resultant",
            (polynomial_key(left), polynomial_key(right)),
        )
    return tuple(claims)


def build_lazard_complete_projection(polys, variables):
    """Build a full Lazard tower, including contents through the factor basis.

    Every level is an irreducible/squarefree factor basis. Factors independent
    of the current variable carry its content downward. This is the full
    leading/trailing-coefficient operator, with no equational-constraint pruning.
    """
    variables = tuple(variables)
    if not variables:
        raise ValueError("at least one variable is required")
    family = _cad_input_polynomials(polys, variables)
    key = ("lazard-complete", polynomial_family_key(family), variables)
    cached = PROJECTION_TOWERS.get(key)
    if cached is not None:
        return cached
    current = squarefree_basis(family)
    levels = {len(variables): ProjectionLevel(len(variables), variables[-1], current)}
    for level in range(len(variables), 1, -1):
        claims = lazard_projection_entries(current, variables[level - 1], variables[: level - 1])
        current = squarefree_basis(claim.poly for claim in claims)
        levels[level - 1] = ProjectionLevel(level - 1, variables[level - 2], current, claims)
    tower = ProjectionTower(
        variables,
        tuple(levels[level] for level in range(1, len(variables) + 1)),
        levels[len(variables)].polynomials,
        {"projection": "lazard-complete", "complete": True, "provenance": True},
    )
    PROJECTION_TOWERS.put(key, tower)
    return tower
