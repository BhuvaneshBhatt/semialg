from __future__ import annotations

import pytest
import sympy as sp

from semialg.algebraic.rational_univariate import (
    compute_rational_univariate_representation,
    solve_zero_dimensional_system_with_rur,
)


@pytest.mark.parametrize(
    ("equations", "quotient_dimension", "solution_count"),
    [
        (lambda x, y: (x**2, y**2), 4, 1),
        (lambda x, y: (x**2, y**3), 6, 1),
        (lambda x, y: ((x - 1) ** 3, y - x), 3, 1),
        (lambda x, y: ((x**2 - 2) ** 3, y - x), 6, 2),
    ],
)
def test_nonreduced_rur_separates_quotient_multiplicity_from_geometric_branches(
    equations, quotient_dimension, solution_count
):
    x, y = sp.symbols("x y", real=True)
    system = equations(x, y)

    representation = compute_rational_univariate_representation(system, (x, y))
    solutions = solve_zero_dimensional_system_with_rur(system, (x, y), real=True)

    assert representation.quotient_dimension == quotient_dimension
    assert representation.geometric_solution_count == solution_count
    assert representation.solution_count == solution_count
    assert len(solutions) == solution_count


@pytest.mark.parametrize("alpha", [sp.sqrt(2), sp.sqrt(3), sp.sqrt(2) + sp.sqrt(3)])
def test_nonreduced_rur_over_algebraic_field_preserves_nilpotent_dimension(alpha):
    x, y = sp.symbols("x y", real=True)
    system = ((x - alpha) ** 2, (y - x) ** 2)

    representation = compute_rational_univariate_representation(system, (x, y))
    solutions = solve_zero_dimensional_system_with_rur(system, (x, y), real=True)

    assert getattr(representation.defining_polynomial.domain, "is_AlgebraicField", False)
    assert representation.quotient_dimension == 4
    assert representation.geometric_solution_count == 1
    assert solutions == ((alpha, alpha),)


def test_nonreduced_rur_cache_survives_rational_algebraic_rational_transition():
    x, y = sp.symbols("x y", real=True)
    rational = ((x - 1) ** 3, y - x)
    algebraic = (sp.sqrt(3) * (x - 1) ** 3, sp.sqrt(3) * (y - x))

    first = compute_rational_univariate_representation(rational, (x, y))
    middle = compute_rational_univariate_representation(algebraic, (x, y))
    final = compute_rational_univariate_representation(rational, (x, y))

    assert str(first.defining_polynomial.domain) == "QQ"
    assert getattr(middle.defining_polynomial.domain, "is_AlgebraicField", False)
    assert str(final.defining_polynomial.domain) == "QQ"
    assert final.quotient_dimension == 3
    assert final.geometric_solution_count == 1


def test_nonradical_two_branch_system_returns_distinct_geometric_solutions_once():
    x, y = sp.symbols("x y", real=True)
    system = ((x**2 - 2) ** 2, y - x)

    solutions = solve_zero_dimensional_system_with_rur(system, (x, y), real=True)

    assert solutions == ((-sp.sqrt(2), -sp.sqrt(2)), (sp.sqrt(2), sp.sqrt(2)))
