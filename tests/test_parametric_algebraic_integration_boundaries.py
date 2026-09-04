import pytest
import sympy as sp

from semialg import integrate_over_region, semialgebraic_measure

x, a = sp.symbols("x a", real=True)


@pytest.fixture(scope="module")
def quadratic_integrals():
    return {
        p: integrate_over_region(x**p, x**2 <= a, [x], parameters=[a], return_stratified=True)
        for p in range(5)
    }


@pytest.mark.parametrize(
    "power,param,expected",
    [
        (0, -1, 0),
        (0, 0, 0),
        (0, 1, 2),
        (0, 4, 4),
        (1, 1, 0),
        (1, 4, 0),
        (2, 1, sp.Rational(2, 3)),
        (2, 4, sp.Rational(16, 3)),
        (3, 1, 0),
        (3, 4, 0),
        (4, 1, sp.Rational(2, 5)),
        (4, 4, sp.Rational(64, 5)),
    ],
)
def test_quadratic_algebraic_endpoints_across_cells(quadratic_integrals, power, param, expected):
    assert sp.simplify(quadratic_integrals[power].select({a: param}) - expected) == 0


@pytest.mark.parametrize("strict", [False, True])
def test_measure_is_unchanged_by_removing_quadratic_boundary(strict):
    relation = x**2 < a if strict else x**2 <= a
    result = semialgebraic_measure(relation, [x], parameters=[a], return_stratified=True)
    assert result.select({a: 0}) == 0
    assert result.select({a: 9}) == 6


@pytest.mark.parametrize(
    "param,expected",
    [(-1, 0), (0, 0), (1, 2), (16, 4)],
)
def test_quartic_root_endpoint_measure(param, expected):
    result = semialgebraic_measure(x**4 <= a, [x], parameters=[a], return_stratified=True)
    assert sp.simplify(result.select({a: param}) - expected) == 0


def test_colliding_roots_at_discriminant_boundary_are_counted_once():
    result = semialgebraic_measure(
        (x**2 - a) ** 2 <= 0, [x], parameters=[a], return_stratified=True
    )
    assert result.select({a: -1}) == 0
    assert result.select({a: 0}) == 0
    assert result.select({a: 4}) == 0
