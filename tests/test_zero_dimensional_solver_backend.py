import pytest
import sympy as sp

from semialg import is_zero_dimensional
from semialg.algebraic.rational_univariate import RationalUnivariateError
from semialg.solve.zero_dimensional import solve_zero_dimensional_system


def test_public_zero_dimensional_solver_uses_rur_and_filters_inequalities():
    x, y = sp.symbols("x y", real=True)

    result = solve_zero_dimensional_system(
        [sp.Eq(x**2 + y**2, 1), sp.Eq(x - y, 0)],
        inequalities=x > 0,
        variables=[x, y],
        return_result=True,
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
        solve_zero_dimensional_system(
            [x**2 + y**2 - 1], variables=[x, y], backend="rur", return_result=True
        )


def test_rur_nonradical_metadata():
    x, y = sp.symbols("x y")

    result = solve_zero_dimensional_system([x**2, y - 1], variables=[x, y], return_result=True)

    assert result.points == ((0, 1),)
    assert result.representation is not None
    assert result.representation.dimension == 2
    assert result.representation.solution_count == 1
    assert result.representation.separating_linear_form == x + y


def test_zero_dimensional_solver_uses_radical_simplification_before_enumeration():
    x = sp.symbols("x", real=True)
    result = solve_zero_dimensional_system(
        [x**2], inequalities=sp.Ne(x, 0), variables=[x], return_result=True
    )
    assert result.points == tuple()
    assert result.status == "unsat"
    assert result.backend == "groebner_ideal+rational_univariate"
    assert result.ideal_analysis is not None
    assert result.ideal_analysis.quotient_dimension == 2
