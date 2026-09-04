from __future__ import annotations

import sympy as sp

from semialg import SemialgebraicRegion, simplify_region


def test_boolean_idempotence_and_absorption_simplify_semantically():
    x = sp.symbols("x", real=True)
    a = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))
    b = SemialgebraicRegion(sp.And(x >= 0, x <= 2), (x,))

    assert simplify_region(a.union(a), exact=True).equals_region(a)
    assert simplify_region(a.intersection(a), exact=True).equals_region(a)
    assert simplify_region(a.union(a.intersection(b)), exact=True).equals_region(a)


def test_exact_union_and_intersection_remove_subset_redundancy():
    x = sp.symbols("x", real=True)
    small = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))
    large = SemialgebraicRegion(sp.And(x >= 0, x <= 2), (x,))

    assert simplify_region(small.union(large), exact=True).equals_region(large)
    assert simplify_region(small.intersection(large), exact=True).equals_region(small)


def test_topological_idempotence_on_interval():
    x = sp.symbols("x", real=True)
    open_interval = SemialgebraicRegion(sp.And(x > 0, x < 1), (x,))
    assert open_interval.interior().interior().equals_region(open_interval.interior())
    assert open_interval.closure().closure().equals_region(open_interval.closure())
