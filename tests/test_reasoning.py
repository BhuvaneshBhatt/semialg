import sympy as sp

from semialg import (
    prove_negative,
    prove_nonnegative,
    prove_nonpositive,
    prove_positive,
    region_bounded,
    region_closed,
    region_compact,
    region_disjoint,
    region_equal,
    region_subset,
    simplify_system,
    simplify_under_assumptions,
)


def test_simplify_system_removes_redundant_constraints():
    x = sp.symbols("x", real=True)
    assert simplify_system([x > 0, x >= 0], [x]) == (x > 0)


def test_simplify_system_detects_contradiction():
    x, y = sp.symbols("x y", real=True)
    assert simplify_system([x**2 < 0, y > 0], [x, y]) == sp.false


def test_inequality_provers():
    x = sp.symbols("x", real=True)
    assert prove_positive(x**2 + 1, [x])
    assert prove_nonnegative((x - 1) ** 2, [x])
    assert prove_negative(-(x**2) - 1, [x])
    assert prove_nonpositive(-(x**2), [x])
    assert prove_nonnegative(x, [x], assumptions=x >= 0)


def test_region_relations():
    x, y = sp.symbols("x y", real=True)
    assert region_subset(x > 1, x > 0, [x])
    assert region_equal(x**2 <= 1, sp.And(x >= -1, x <= 1), [x])
    assert region_disjoint(x < 0, x > 0, [x])
    assert region_subset(x**2 + y**2 <= 1, x**2 + y**2 <= 4, [x, y])


def test_boundedness_closedness_compactness():
    x, y = sp.symbols("x y", real=True)
    assert region_bounded(sp.And(x >= 0, x <= 1), [x])
    assert not region_bounded(x >= 0, [x])
    assert region_bounded(x**2 + y**2 <= 1, [x, y])
    assert region_closed(sp.And(x >= 0, x <= 1), [x])
    assert not region_closed(sp.And(x > 0, x < 1), [x])
    assert region_compact(sp.And(x >= 0, x <= 1), [x])
    assert not region_compact(sp.And(x > 0, x < 1), [x])


def test_simplify_under_assumptions_abs_sqrt_minmax_piecewise():
    x, y = sp.symbols("x y", real=True)
    assert simplify_under_assumptions(sp.Abs(x), x >= 0, [x]) == x
    assert simplify_under_assumptions(sp.Abs(x), x <= 0, [x]) == -x
    assert simplify_under_assumptions(sp.sqrt(x**2), x >= 0, [x]) == x
    assert simplify_under_assumptions(sp.sqrt((x - 1) ** 2), x >= 1, [x]) == x - 1
    assert simplify_under_assumptions(sp.Max(x, y), x >= y, [x, y]) == x
    assert simplify_under_assumptions(sp.Min(x, y), x >= y, [x, y]) == y
    expr = sp.Piecewise((x, x >= 0), (-x, True))
    assert simplify_under_assumptions(expr, x >= 0, [x]) == x
