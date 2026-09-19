"""Behavioral calls through the package root for central user workflows."""

import sympy as sp

import semialg


def test_root_decision_and_solving_workflow() -> None:
    x = sp.Symbol("x", real=True)
    assert semialg.is_satisfiable(x**2 < 0, [x]) is False
    assert semialg.equivalent(x**2 <= 0, sp.Eq(x, 0), [x]) is True
    solution = semialg.solve_semialgebraic(sp.And(x >= 0, x <= 1), [x])
    assert solution.satisfiable is True
    assert solution.dimension == 1
    assert solution.samples == ({x: sp.Rational(1, 2)},)
    assert semialg.find_instance(sp.And(x > 0, x < 2), [x]) == {x: 1}


def test_root_function_analysis_workflow() -> None:
    x, t = sp.symbols("x t", real=True)
    assert semialg.function_domain(sp.sqrt(x), [x]) == (x >= 0)
    assert semialg.equivalent(
        semialg.function_range(x**2, variables=[x], value_symbol=t), t >= 0, [t]
    )
    assert semialg.is_real_valued(sp.sqrt(x), [x], assumptions=x >= 0) is True


def test_root_region_measure_and_integration_workflow() -> None:
    x, y = sp.symbols("x y", real=True)
    region = sp.And(x >= 0, x <= 1, y >= 0, y <= x)
    assert semialg.integrate_over_region(1, region, [x, y]) == sp.Rational(1, 2)
    assert semialg.region_dimension(region, [x, y]) == 2
    assert semialg.is_bounded(region, [x, y]) is True
