from __future__ import annotations

import pytest
import sympy as sp

from semialg import is_convex
from semialg.geometry_queries import (
    _is_affine_polyhedron,
    _is_convex_by_definition,
    _is_convex_polynomial_intersection,
)
from semialg.normalization import normalize_formula


def _formula(expr):
    return normalize_formula(expr)


def test_affine_polyhedron_fast_path_matches_definition():
    x, y = sp.symbols("x y", real=True)
    region = sp.And(x >= 0, y >= 0, x + y <= 1)
    formula = _formula(region)

    assert _is_affine_polyhedron(formula, (x, y))
    assert is_convex(region, (x, y))


def test_open_affine_polyhedron_is_convex():
    x, y = sp.symbols("x y", real=True)
    region = sp.And(x > 0, y > 0, x + y < 1)
    formula = _formula(region)

    assert _is_affine_polyhedron(formula, (x, y))
    assert is_convex(region, (x, y))


def test_lower_dimensional_affine_polyhedron_is_convex():
    x, y = sp.symbols("x y", real=True)
    region = sp.And(sp.Eq(y, 0), x >= -1, x <= 2)
    formula = _formula(region)

    assert _is_affine_polyhedron(formula, (x, y))
    assert is_convex(region, (x, y))


@pytest.mark.parametrize(
    "region",
    [
        lambda x, y: x**4 + y**4 <= 1,
        lambda x, y: (x**2 + y**2) ** 2 <= 1,
        lambda x, y: sp.And(x**4 + y**4 <= 1, x + y >= -1),
        lambda x, y: -(x**4 + y**4) >= -1,
        lambda x, y: sp.And(x**4 + y**4 < 1, x - y <= sp.Rational(1, 2)),
    ],
)
def test_convex_polynomial_sublevel_fast_path_matches_definition(region):
    x, y = sp.symbols("x y", real=True)
    condition = region(x, y)
    formula = _formula(condition)

    assert _is_convex_polynomial_intersection(formula, (x, y))
    assert is_convex(condition, (x, y))


def test_convex_polynomial_sublevel_with_affine_equality():
    x, y = sp.symbols("x y", real=True)
    region = sp.And(x**4 + y**4 <= 1, sp.Eq(x + y, 0))
    formula = _formula(region)

    assert _is_convex_polynomial_intersection(formula, (x, y))
    assert is_convex(region, (x, y))


def test_nonconvex_polynomial_superlevel_is_not_certified_by_hessian_path():
    x, y = sp.symbols("x y", real=True)
    region = x**2 + y**2 >= 1
    formula = _formula(region)

    assert not _is_convex_polynomial_intersection(formula, (x, y))
    assert bool(region.subs({x: -1, y: 0}))
    assert bool(region.subs({x: 1, y: 0}))
    assert not bool(region.subs({x: 0, y: 0}))


def test_disconnected_union_is_nonconvex_by_complete_fallback():
    x, y = sp.symbols("x y", real=True)
    left = sp.And(x >= -2, x <= -1, y >= -1, y <= 1)
    right = sp.And(x >= 1, x <= 2, y >= -1, y <= 1)
    region = sp.Or(left, right)
    formula = _formula(region)

    assert not _is_affine_polyhedron(formula, (x, y))
    assert not _is_convex_polynomial_intersection(formula, (x, y))
    assert bool(region.subs({x: -sp.Rational(3, 2), y: 0}))
    assert bool(region.subs({x: sp.Rational(3, 2), y: 0}))
    assert not bool(region.subs({x: 0, y: 0}))


def test_punctured_disk_is_nonconvex():
    x, y = sp.symbols("x y", real=True)
    region = sp.And(x**2 + y**2 <= 1, sp.Ne(x**2 + y**2, 0))
    formula = _formula(region)

    assert not _is_convex_polynomial_intersection(formula, (x, y))
    assert bool(region.subs({x: -1, y: 0}))
    assert bool(region.subs({x: 1, y: 0}))
    assert not bool(region.subs({x: 0, y: 0}))


def test_nonlinear_singleton_falls_back_and_is_convex():
    x, y = sp.symbols("x y", real=True)
    region = sp.Eq(x**2 + y**2, 0)
    formula = _formula(region)

    assert not _is_affine_polyhedron(formula, (x, y))
    assert not _is_convex_polynomial_intersection(formula, (x, y))
    assert is_convex(region, (x, y))


def test_affine_change_of_coordinates_preserves_convexity():
    x, y = sp.symbols("x y", real=True)
    u = 2 * x + y
    v = x - y
    transformed_disk = u**2 + v**2 <= 1

    assert is_convex(transformed_disk, (x, y))


@pytest.mark.slow
def test_quantified_definition_accepts_interval():
    x = sp.symbols("x", real=True)
    formula = _formula(sp.And(x >= -1, x <= 2))
    assert _is_convex_by_definition(formula, (x,))


@pytest.mark.slow
def test_quantified_definition_rejects_two_point_set():
    x = sp.symbols("x", real=True)
    formula = _formula(sp.Or(sp.Eq(x, -1), sp.Eq(x, 1)))
    assert _is_convex_by_definition(formula, (x,)) is False
