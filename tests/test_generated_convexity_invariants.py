from __future__ import annotations

import pytest
import sympy as sp

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg import is_convex
from semialg.geometry_queries import (
    _is_convex_polynomial_intersection,
)
from semialg.normalization import normalize_formula

small = st.integers(min_value=-2, max_value=2)
positive = st.integers(min_value=1, max_value=3)


@settings(max_examples=20, deadline=None)
@given(a=positive, b=positive, c=small, d=small)
def test_generated_separable_quartic_sublevels_match_definition(a, b, c, d):
    x, y = sp.symbols("x y", real=True)
    poly = a * (x - c) ** 4 + b * (y - d) ** 4
    region = poly <= 4
    formula = normalize_formula(region)

    assert _is_convex_polynomial_intersection(formula, (x, y))
    assert is_convex(region, (x, y))


@settings(max_examples=20, deadline=None)
@given(
    a=positive,
    b=positive,
    c=small,
    d=small,
    rhs=st.integers(min_value=0, max_value=4),
)
def test_generated_convex_quadratic_polyhedra_match_definition(a, b, c, d, rhs):
    x, y = sp.symbols("x y", real=True)
    region = sp.And(a * x**2 + b * y**2 <= rhs + 1, c * x + d * y <= rhs + 2)
    formula = normalize_formula(region)

    assert _is_convex_polynomial_intersection(formula, (x, y))
    assert is_convex(region, (x, y))


@settings(max_examples=12, deadline=None)
@given(offset=st.integers(min_value=1, max_value=4))
def test_generated_separated_box_unions_are_nonconvex(offset):
    x, y = sp.symbols("x y", real=True)
    left = sp.And(x >= -offset - 2, x <= -offset - 1, y >= -1, y <= 1)
    right = sp.And(x >= offset + 1, x <= offset + 2, y >= -1, y <= 1)
    region = sp.Or(left, right)

    assert is_convex(region, (x, y)) is False


@settings(max_examples=16, deadline=None)
@given(shift=small, scale=positive)
def test_generated_one_dimensional_interval_rewrites_remain_convex(shift, scale):
    x = sp.symbols("x", real=True)
    left = scale * (x - shift) >= -scale
    right = scale * (x - shift) <= 2 * scale
    region = sp.And(left, right)

    assert is_convex(region, (x,)) is True


@settings(max_examples=16, deadline=None)
@given(offset=st.integers(min_value=1, max_value=4))
def test_generated_one_dimensional_two_ray_sets_are_nonconvex(offset):
    x = sp.symbols("x", real=True)
    region = sp.Or(x <= -offset, x >= offset)

    assert is_convex(region, (x,)) is False


@settings(max_examples=12, deadline=None)
@given(bound=st.integers(min_value=2, max_value=6))
def test_generated_domain_relative_quartics_are_certified(bound):
    x, y = sp.symbols("x y", real=True)
    expr = x**4 - x**2 + y**2
    region = sp.And(x >= 1, expr <= bound)

    assert is_convex(region, (x, y)) is True
