import pytest
import sympy as sp

from semialg import is_zero_dimensional, solve_zero_dimensional_system
from semialg.algebraic.rational_univariate import RationalUnivariateError


def test_public_zero_dimensional_solver_uses_rur_and_filters_inequalities():
    x, y = sp.symbols("x y", real=True)

    result = solve_zero_dimensional_system(
        [sp.Eq(x**2 + y**2, 1), sp.Eq(x - y, 0)],
        inequalities=x > 0,
        vars=[x, y],
    )

    assert result.backend == "rational_univariate"
    assert result.status == "satisfied"
    assert len(result.points) == 1
    assert sp.simplify(result.points[0][0] - sp.sqrt(2) / 2) == 0
    assert sp.simplify(result.points[0][1] - sp.sqrt(2) / 2) == 0
    assert result.representation is not None
    assert result.assignments[0][x] == result.points[0][0]


def test_public_zero_dimensional_detection_rejects_positive_dimensional_curve():
    x, y = sp.symbols("x y")

    assert is_zero_dimensional([x**2 + y**2 - 1, x - y], [x, y])
    assert not is_zero_dimensional([x**2 + y**2 - 1], [x, y])

    with pytest.raises(RationalUnivariateError, match="zero-dimensional"):
        solve_zero_dimensional_system([x**2 + y**2 - 1], vars=[x, y], backend="rur")


def test_rur_backend_exposes_quotient_and_geometric_solution_metadata_for_nonradical_system():
    x, y = sp.symbols("x y")

    result = solve_zero_dimensional_system([x**2, y - 1], vars=[x, y])

    assert result.points == ((0, 1),)
    assert result.representation is not None
    assert result.representation.dimension == 2
    assert result.representation.solution_count == 1
    assert result.representation.separating_linear_form == x + y
