from __future__ import annotations

import sympy as sp

from semialg import SemialgebraicRegion


def test_measure_is_integral_of_one():
    x, y = sp.symbols("x y", real=True)
    region = SemialgebraicRegion(sp.And(x >= 0, x <= 2, y >= -1, y <= 1), (x, y))
    assert sp.simplify(region.measure() - region.integrate(1)) == 0


def test_measure_inclusion_exclusion_for_overlapping_intervals():
    x = sp.symbols("x", real=True)
    a = SemialgebraicRegion(sp.And(x >= 0, x <= 2), (x,))
    b = SemialgebraicRegion(sp.And(x >= 1, x <= 3), (x,))
    union = a.union(b)
    inter = a.intersection(b)
    assert sp.simplify(union.measure() - (a.measure() + b.measure() - inter.measure())) == 0


def test_cad_region_measure_reuses_existing_decomposition():
    x = sp.symbols("x", real=True)
    region = SemialgebraicRegion(sp.And(x >= -2, x <= 3), (x,))
    cad = region.as_cad_region()
    result = cad.result
    assert sp.simplify(cad.measure() - 5) == 0
    assert cad.result is result
    assert region.ensure_cad() is result
