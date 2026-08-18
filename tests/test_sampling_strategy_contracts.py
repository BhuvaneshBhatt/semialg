from __future__ import annotations

import sympy as sp

from semialg import sample_point, sample_points


def test_sample_points_grid_strategy_uses_bounds_and_rationals() -> None:
    x, y = sp.symbols("x y", real=True)
    pts = sample_points(
        (x >= 0) & (y >= 0) & (x + y <= 1),
        [x, y],
        count=3,
        strategy="grid",
        bounds=[(0, 1), (0, 1)],
        grid_resolution=3,
    )
    assert len(pts) == 3
    assert all(point[x].is_Rational and point[y].is_Rational for point in pts)
    assert all(bool(((x >= 0) & (y >= 0) & (x + y <= 1)).subs(point)) for point in pts)


def test_sample_points_random_strategy_is_seeded() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = x**2 + y**2 < 1
    kwargs = dict(
        variables=[x, y],
        count=2,
        strategy="random",
        bounds=[(-1, 1), (-1, 1)],
        seed=123,
        exact=False,
        random_attempts=200,
    )
    pts1 = sample_points(formula, **kwargs)
    pts2 = sample_points(formula, **kwargs)
    assert pts1 == pts2
    assert len(pts1) == 2
    assert all(point[x].is_Float and point[y].is_Float for point in pts1)
    assert all(bool(formula.subs(point)) for point in pts1)


def test_sample_point_representative_finds_exact_irrational_rur_witness() -> None:
    x = sp.symbols("x", real=True)
    point = sample_point(sp.Eq(x**2, 2), [x], strategy="representative")
    assert point is not None
    assert sp.simplify(point[x] ** 2 - 2) == 0


def test_sample_points_rational_strategy_does_not_use_float_randomness() -> None:
    x = sp.symbols("x", real=True)
    pts = sample_points(x**2 > 1, [x], count=2, strategy="rational")
    assert len(pts) == 2
    assert all(point[x].is_Rational for point in pts)
