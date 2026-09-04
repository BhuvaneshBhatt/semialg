import sympy as sp

from semialg import integrate_over_region, semialgebraic_measure


def test_parametric_quadratic_root_endpoints_integrate_exactly():
    x, a = sp.symbols("x a", real=True)
    result = integrate_over_region(
        1,
        x**2 <= a,
        [x],
        parameters=[a],
        return_stratified=True,
    )
    assert result.method == "parametric_algebraic_root_cad_integration"
    assert result.select({a: -1}) == 0
    assert result.select({a: 0}) == 0
    assert result.select({a: 4}) == 4
    assert result.select({a: 9}) == 6


def test_parametric_quadratic_measure_uses_algebraic_root_endpoints():
    x, a = sp.symbols("x a", real=True)
    result = semialgebraic_measure(
        x**2 < a,
        [x],
        parameters=[a],
        return_stratified=True,
    )
    assert result.select({a: 1}) == 2
    assert result.select({a: 16}) == 8
