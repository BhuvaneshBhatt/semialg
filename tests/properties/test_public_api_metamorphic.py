"""Generated transformations with explicit exact mathematical expectations."""

from __future__ import annotations

import pytest
import sympy as sp
from hypothesis import example, given, settings
from hypothesis import strategies as st

from semialg import (
    function_range,
    integrate_over_region,
    region_complement,
    region_difference,
    region_intersection,
    region_union,
    scale,
    semialgebraic_image,
    semialgebraic_maximize,
    semialgebraic_minimize,
    semialgebraic_preimage,
    translate,
)

pytestmark = pytest.mark.slow

X, U = sp.symbols("x u", real=True)
SMALL = st.integers(-3, 3)
NONZERO = st.sampled_from((-3, -2, -1, 1, 2, 3))


@st.composite
def intervals(draw):
    left = draw(SMALL)
    width = draw(st.integers(0, 4))
    return sp.Interval(left, left + width, draw(st.booleans()), draw(st.booleans()))


def _real_set(formula):
    return formula.as_set() & sp.S.Reals


@settings(max_examples=24, deadline=None)
@given(left=intervals(), right=intervals())
@example(left=sp.Interval.open(0, 1), right=sp.Interval(1, 2))
@example(left=sp.S.EmptySet, right=sp.FiniteSet(0))
def test_boolean_region_algebra_matches_exact_sets(left, right):
    first, second = left.as_relational(X), right.as_relational(X)
    assert _real_set(region_union(first, second)) == left | right
    assert _real_set(region_intersection(first, second)) == left & right
    assert _real_set(region_difference(first, second)) == left - right
    complement = region_complement(first)
    assert _real_set(complement) == sp.S.Reals - left
    assert _real_set(region_complement(complement)) == left
    de_morgan = region_intersection(complement, region_complement(second))
    assert _real_set(de_morgan) == sp.S.Reals - (left | right)


@settings(max_examples=18, deadline=None)
@given(
    left=SMALL,
    width=st.integers(1, 4),
    factor=NONZERO,
    offset=SMALL,
    left_open=st.booleans(),
    right_open=st.booleans(),
)
@example(left=0, width=2, factor=-2, offset=1, left_open=True, right_open=False)
def test_affine_image_preimage_and_inverse_preserve_endpoints(
    left, width, factor, offset, left_open, right_open
):
    source = sp.Interval(left, left + width, left_open, right_open)
    domain = source.as_relational(X)
    ends = (factor * left + offset, factor * (left + width) + offset)
    expected = (
        sp.Interval(*ends, left_open, right_open)
        if factor > 0
        else sp.Interval(ends[1], ends[0], right_open, left_open)
    )
    image = semialgebraic_image((factor * X + offset,), domain, (X,), image_variables=(U,))
    assert _real_set(image) == expected
    composed = translate(scale(domain, factor, (X,)), (offset,), (X,))
    assert _real_set(composed) == expected
    preimage = semialgebraic_preimage((factor * X + offset,), image, (X,), target_variables=(U,))
    assert _real_set(preimage) == source
    inverse = scale(translate(composed, (-offset,), (X,)), sp.Rational(1, factor), (X,))
    assert _real_set(inverse) == source


@settings(max_examples=18, deadline=None)
@given(left=SMALL, width=st.integers(1, 4), factor=NONZERO, offset=SMALL, degree=st.integers(0, 3))
@example(left=-1, width=3, factor=-2, offset=1, degree=3)
def test_integral_affine_change_of_variables(left, width, factor, offset, degree):
    domain = sp.Interval(left, left + width).as_relational(X)
    low, high = sorted((factor * left + offset, factor * (left + width) + offset))
    target = sp.Interval(low, high).as_relational(U)
    expected = (sp.Rational(high) ** (degree + 1) - sp.Rational(low) ** (degree + 1)) / (degree + 1)
    direct = integrate_over_region(U**degree, target, (U,))
    pulled_back = integrate_over_region(abs(factor) * (factor * X + offset) ** degree, domain, (X,))
    assert direct == pulled_back == expected


@settings(max_examples=20, deadline=None)
@given(left=SMALL, width=st.integers(2, 5), degree=st.integers(0, 3), coefficient=SMALL)
def test_integral_additivity_and_linearity(left, width, degree, coefficient):
    split = sp.Rational(2 * left + width, 2)
    whole = sp.Interval(left, left + width).as_relational(X)
    lower = sp.Interval(left, split).as_relational(X)
    upper = sp.Interval(split, left + width).as_relational(X)
    moment = (sp.Rational(left + width) ** (degree + 1) - sp.Rational(left) ** (degree + 1)) / (
        degree + 1
    )
    total = integrate_over_region(coefficient * X**degree + 1, whole, (X,))
    parts = integrate_over_region(coefficient * X**degree + 1, lower, (X,)) + integrate_over_region(
        coefficient * X**degree + 1, upper, (X,)
    )
    assert total == parts == coefficient * moment + width


@settings(max_examples=18, deadline=None)
@given(center=SMALL, radius=st.integers(1, 3), factor=NONZERO, offset=SMALL)
def test_objective_scaling_swaps_extrema_and_preserves_optimizers(center, radius, factor, offset):
    domain = sp.Interval(center - radius, center + radius).as_relational(X)
    objective = factor * (X - center) ** 2 + offset
    low, high = sorted((offset, factor * radius**2 + offset))
    minimum = semialgebraic_minimize(objective, domain, (X,), return_result=True)
    maximum = semialgebraic_maximize(objective, domain, (X,), return_result=True)
    image = function_range(objective, domain, (X,), return_result=True)
    assert minimum.value == low and maximum.value == high
    assert minimum.attained and maximum.attained
    assert minimum.certified and maximum.certified
    assert _real_set(image.formula) == sp.Interval(low, high)
    for result in (minimum, maximum):
        assert result.points
        for point in result.points:
            assert domain.subs(point) is sp.true
            assert objective.subs(point) == result.value
    center_result = minimum if factor > 0 else maximum
    assert {point[X] for point in center_result.points} == {center}
