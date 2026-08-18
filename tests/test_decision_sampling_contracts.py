from __future__ import annotations

import sympy as sp

from semialg import sample_point, sample_points, solve_semialgebraic


def test_sample_points_rational_strategy_returns_exact_valid_points():
    x = sp.symbols("x", real=True)

    points = sample_points((x >= 0) & (x <= 1), [x], count=3, strategy="rational")

    assert len(points) == 3
    assert all(point[x].is_Rational for point in points)
    assert all(bool(((x >= 0) & (x <= 1)).subs(point)) for point in points)


def test_sample_points_grid_strategy_uses_bounds():
    x, y = sp.symbols("x y", real=True)
    formula = (x >= 1) & (x <= 2) & (y >= -1) & (y <= 0)

    points = sample_points(
        formula,
        [x, y],
        count=2,
        strategy="grid",
        bounds=[(1, 2), (-1, 0)],
        grid_resolution=3,
    )

    assert len(points) == 2
    assert all(1 <= point[x] <= 2 and -1 <= point[y] <= 0 for point in points)


def test_sample_points_random_strategy_is_seeded_and_can_return_floats():
    x = sp.symbols("x", real=True)
    formula = (x >= 0) & (x <= 1)

    left = sample_points(formula, [x], count=3, strategy="random", exact=False, seed=123)
    right = sample_points(formula, [x], count=3, strategy="random", exact=False, seed=123)

    assert left == right
    assert len(left) == 3
    assert all(isinstance(point[x], sp.Float) for point in left)


def test_sample_point_representative_still_uses_rur_for_irrational_finite_witness():
    x = sp.symbols("x", real=True)

    point = sample_point(sp.Eq(x**2, 2), [x], strategy="representative")

    assert point is not None
    assert sp.simplify(point[x] ** 2 - 2) == 0


def test_solve_semialgebraic_accepts_random_sampling_strategy():
    x = sp.symbols("x", real=True)

    result = solve_semialgebraic((x >= 0) & (x <= 1), [x], count=2, strategy="random")

    assert result.satisfiable
    assert len(result.samples) == 2
    assert all(bool(((x >= 0) & (x <= 1)).subs(point)) for point in result.samples)
