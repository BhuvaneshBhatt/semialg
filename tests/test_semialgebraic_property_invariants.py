import random

import pytest
import sympy as sp

from semialg import (
    affine_transform,
    centroid,
    coordinate_range,
    is_closed,
    is_compact,
    is_disjoint,
    is_empty,
    is_equal,
    is_open,
    is_subset,
    region_boundary,
    region_closure,
    region_dimension,
    region_interior,
    semialgebraic_measure,
    translate,
)

x, y = sp.symbols("x y", real=True)


def _interval(lo, hi):
    return sp.And(x >= lo, x <= hi)


RNG = random.Random(20260827)
INTERVALS = []
for _ in range(12):
    lo = RNG.randint(-8, 3)
    hi = RNG.randint(lo + 1, lo + 8)
    INTERVALS.append((sp.Integer(lo), sp.Integer(hi)))


@pytest.mark.parametrize("lo,hi", INTERVALS)
def test_generated_interval_measure_centroid_range_and_compactness(lo, hi):
    region = _interval(lo, hi)
    assert semialgebraic_measure(region, [x]) == hi - lo
    assert centroid(region, [x])[x] == sp.Rational(lo + hi, 2)
    assert is_compact(region, [x])
    rng = coordinate_range(region, x, [x])
    value_symbol = next(sym for sym in rng.free_symbols if sym != x)
    assert sp.simplify(rng.subs(value_symbol, lo)) is sp.true
    assert sp.simplify(rng.subs(value_symbol, hi)) is sp.true
    assert sp.simplify(rng.subs(value_symbol, hi + 1)) is sp.false


@pytest.mark.parametrize("lo,hi", INTERVALS[:8])
def test_generated_set_algebra_identities(lo, hi):
    a = _interval(lo, hi)
    b = _interval(lo - 1, hi + 1)
    assert is_equal(sp.Or(a, a), a, [x])
    assert is_equal(sp.And(a, a), a, [x])
    assert is_subset(a, b, [x])
    assert is_empty(sp.And(a, sp.Not(a)), [x])
    assert is_disjoint(a, x > hi, [x])


@pytest.mark.parametrize("lo,hi", INTERVALS[:8])
def test_generated_topology_identities(lo, hi):
    closed = _interval(lo, hi)
    opened = sp.And(x > lo, x < hi)
    assert is_equal(
        region_closure(region_closure(opened, [x]), [x]), region_closure(opened, [x]), [x]
    )
    assert is_equal(
        region_interior(region_interior(closed, [x]), [x]), region_interior(closed, [x]), [x]
    )
    expected_boundary = sp.And(region_closure(closed, [x]), sp.Not(region_interior(closed, [x])))
    assert is_equal(region_boundary(closed, [x]), expected_boundary, [x])
    assert is_closed(closed, [x])
    assert is_open(opened, [x])


@pytest.mark.parametrize("lo,hi", INTERVALS[:6])
def test_generated_affine_translation_invariants(lo, hi):
    region = _interval(lo, hi)
    shift = sp.Integer(3)
    moved = translate(region, (shift,), [x])
    assert semialgebraic_measure(moved, [x]) == semialgebraic_measure(region, [x])
    assert region_dimension(moved, [x]) == region_dimension(region, [x])
    assert centroid(moved, [x])[x] == centroid(region, [x])[x] + shift


@pytest.mark.parametrize(
    "bounds",
    [(-2, 1, -3, 2), (0, 2, 1, 4), (-4, -1, -2, 3), (-1, 3, -1, 3)],
)
def test_generated_boxes_have_known_volume_centroid_and_dimension(bounds):
    x0, x1, y0, y1 = map(sp.Integer, bounds)
    box = sp.And(x >= x0, x <= x1, y >= y0, y <= y1)
    assert semialgebraic_measure(box, [x, y]) == (x1 - x0) * (y1 - y0)
    assert centroid(box, [x, y]) == {
        x: sp.Rational(x0 + x1, 2),
        y: sp.Rational(y0 + y1, 2),
    }
    assert region_dimension(box, [x, y]) == 2
    assert is_compact(box, [x, y])


@pytest.mark.slow
def test_nonsingular_affine_image_preserves_full_dimension_of_box():
    box = sp.And(x >= -1, x <= 1, y >= -2, y <= 2)
    image = affine_transform(box, [[1, 1], [0, 2]], [3, -1], [x, y])
    assert region_dimension(image, [x, y]) == 2
    assert is_compact(image, [x, y])
