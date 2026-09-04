import pytest
import sympy as sp

from semialg import (
    affine_transform,
    argmin_set,
    centroid,
    coordinate_range,
    covariance_matrix,
    diameter,
    inertia_tensor,
    intersects,
    is_compact,
    is_disjoint,
    is_equal,
    is_subset,
    minkowski_sum,
    moment_matrix,
    scale,
    support_function,
    translate,
    width,
)

x, y = sp.symbols("x y", real=True)


@pytest.mark.parametrize(
    "lo,hi,shift",
    [
        (-2, 3, 5),
        pytest.param(0, 1, -3, marks=pytest.mark.slow),
        pytest.param(-5, -1, 2, marks=pytest.mark.slow),
        pytest.param(2, 7, -4, marks=pytest.mark.slow),
    ],
)
def test_translation_preserves_diameter_and_shifts_centroid(lo, hi, shift):
    region = sp.And(x >= lo, x <= hi)
    moved = translate(region, (shift,), [x])
    assert diameter(moved, [x]) == diameter(region, [x])
    assert centroid(moved, [x])[x] == centroid(region, [x])[x] + shift


@pytest.mark.parametrize(
    "lo,hi,shift",
    [
        (-2, 3, 5),
        pytest.param(0, 1, -3, marks=pytest.mark.slow),
        pytest.param(-5, -1, 2, marks=pytest.mark.slow),
    ],
)
def test_support_function_translation_identity(lo, hi, shift):
    region = sp.And(x >= lo, x <= hi)
    moved = translate(region, (shift,), [x])
    assert (
        sp.simplify(
            support_function(moved, (1,), [x]) - support_function(region, (1,), [x]) - shift
        )
        == 0
    )


@pytest.mark.parametrize(
    "lo,hi",
    [
        (-2, 3),
        pytest.param(0, 1, marks=pytest.mark.slow),
        pytest.param(-5, -1, marks=pytest.mark.slow),
        pytest.param(2, 7, marks=pytest.mark.slow),
    ],
)
def test_width_equals_two_sided_support_sum(lo, hi):
    region = sp.And(x >= lo, x <= hi)
    expected = support_function(region, (1,), [x]) + support_function(region, (-1,), [x])
    assert sp.simplify(width(region, (1,), [x]) - expected) == 0


@pytest.mark.parametrize(
    "lam",
    [
        2,
        pytest.param(sp.Rational(1, 2), marks=pytest.mark.slow),
        pytest.param(3, marks=pytest.mark.slow),
    ],
)
def test_positive_scaling_scales_diameter_and_width(lam):
    region = sp.And(x >= -1, x <= 2)
    scaled = scale(region, lam, [x])
    assert sp.simplify(diameter(scaled, [x]) - lam * diameter(region, [x])) == 0
    assert sp.simplify(width(scaled, (1,), [x]) - lam * width(region, (1,), [x])) == 0


def test_minkowski_support_function_identity_for_intervals():
    left = sp.And(x >= -1, x <= 2)
    right = sp.And(x >= 3, x <= 5)
    summed = minkowski_sum(left, right, [x])
    assert support_function(summed, (1,), [x]) == support_function(
        left, (1,), [x]
    ) + support_function(right, (1,), [x])


def test_affine_bijection_preserves_subset_relations():
    inner = sp.And(x >= 0, x <= 1)
    outer = sp.And(x >= -1, x <= 2)
    tin = affine_transform(inner, [[3]], [4], [x])
    tout = affine_transform(outer, [[3]], [4], [x])
    assert is_subset(tin, tout, [x])


def test_disjointness_and_intersection_are_complements():
    left = x <= 0
    right = x >= 1
    assert is_disjoint(left, right, [x])
    assert not intersects(left, right, [x])


def test_argmin_set_is_invariant_under_positive_affine_objective_change():
    region = sp.And(x >= -2, x <= 3)
    base = argmin_set((x - 1) ** 2, region, [x])
    transformed = argmin_set(7 * (x - 1) ** 2 + 11, region, [x])
    assert is_equal(base, transformed, [x])


def test_moment_covariance_identity_on_centered_square():
    square = sp.And(x >= -1, x <= 1, y >= -1, y <= 1)
    assert moment_matrix(square, [x, y]) == covariance_matrix(square, [x, y])
    assert inertia_tensor(square, [x, y]) == sp.diag(sp.Rational(4, 3), sp.Rational(4, 3))


@pytest.mark.parametrize("point,expected", [(-2, True), (0, True), (2, True), (3, False)])
def test_coordinate_range_agrees_with_membership(point, expected):
    region = sp.And(x >= -2, x <= 2)
    rng = coordinate_range(region, x, [x])
    value_symbol = next(sym for sym in rng.free_symbols if sym != x)
    assert bool(sp.simplify(rng.subs(value_symbol, point))) is expected


@pytest.mark.parametrize(
    "left,right,subset,equal,disjoint",
    [
        (x <= 0, x <= 1, True, False, False),
        (x < 0, x <= 0, True, False, False),
        (sp.Eq(x, 0), x <= 0, True, False, False),
        (sp.And(x >= 0, x <= 1), sp.And(x >= 0, x <= 1), True, True, False),
        (x < 0, x > 0, False, False, True),
        (x <= -1, x >= 1, False, False, True),
        (sp.false, x >= 0, True, False, True),
        (sp.false, sp.false, True, True, True),
    ],
)
def test_set_relation_algebra_identity_cases(left, right, subset, equal, disjoint):
    assert is_subset(left, right, [x]) is subset
    assert is_equal(left, right, [x]) is equal
    assert is_disjoint(left, right, [x]) is disjoint


def test_compact_interval_remains_compact_under_affine_bijection():
    region = sp.And(x >= -2, x <= 2)
    transformed = affine_transform(region, [[-2]], [5], [x])
    assert is_compact(region, [x])
    assert is_compact(transformed, [x])
