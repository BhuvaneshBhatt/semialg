from __future__ import annotations

import sympy as sp

from semialg import integrate_over_region


def test_zero_dimensional_integral_respects_explicit_bounds():
    x = sp.symbols("x", real=True)
    points = sp.Eq(x**2 - 1, 0)

    assert integrate_over_region(1, points, (x,), bounds={x: (0, 2)}, measure_dimension=0) == 1
    assert integrate_over_region(x, points, (x,), bounds={x: (0, 2)}, measure_dimension=0) == 1


def test_intrinsic_graph_integral_respects_explicit_bounds():
    x, y = sp.symbols("x y", real=True)
    line = sp.Eq(y, x)

    value = integrate_over_region(
        1,
        line,
        (x, y),
        bounds={x: (0, 1)},
        measure_dimension=1,
    )
    assert sp.simplify(value - sp.sqrt(2)) == 0


def test_explicitly_empty_intrinsic_region_has_zero_measure():
    x = sp.symbols("x", real=True)

    assert integrate_over_region(1, sp.false, (x,), measure_dimension="intrinsic") == 0
    assert integrate_over_region(7, sp.false, (x,), measure_dimension="top") == 0


def test_conflicting_circle_constraints_have_zero_intrinsic_measure():
    x, y = sp.symbols("x y", real=True)
    empty = sp.And(sp.Eq(x**2 + y**2, 1), sp.Eq(x**2 + y**2, 4))

    assert integrate_over_region(1, empty, (x, y), measure_dimension=1) == 0
