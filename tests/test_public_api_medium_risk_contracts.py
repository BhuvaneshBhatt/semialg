"""Independent contracts for thin medium-risk root API conveniences."""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib
import sympy as sp

from semialg import (
    Ball,
    Box,
    affine_relative_interior_formula,
    discretize_region_geometry,
    function_convex_partition,
    function_monotonic_partition,
    function_smoothness,
    is_bijective,
    is_function_continuous,
    is_function_smooth,
    is_injective,
    is_surjective,
    is_zero_dimensional,
    matrix_pd_on,
    matrix_psd_on,
    matrix_rank_on,
    matrix_rank_stratification,
    plot_region_geometry,
    plot_solution,
    random_point,
    random_points,
    strict_feasible,
)

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt


def test_function_partition_and_smoothness_conveniences_report_structure() -> None:
    x = sp.Symbol("x", real=True)
    monotonic = function_monotonic_partition(x**2, x)
    convex = function_convex_partition(x**4 - x**2, x)
    smoothness = function_smoothness(sp.Abs(x), x)

    assert tuple(value for value, _ in monotonic) == (
        "strictly_decreasing",
        "strictly_increasing",
    )
    assert tuple(value for value, _ in convex) == ("convex", "concave", "convex")
    assert smoothness.continuous is True
    assert smoothness.smooth is False
    assert smoothness.derivative_exceptions[1] == sp.Eq(x, 0)


def test_mapping_predicates_agree_for_identity_and_square() -> None:
    x = sp.Symbol("x", real=True)

    assert is_injective(x, x) is True
    assert is_surjective(x, x) is True
    assert is_bijective(x, x) is True
    assert is_injective(x**2, x) is False
    assert is_surjective(x**2, x) is False
    assert is_bijective(x**2, x) is False
    assert is_function_continuous(sp.Abs(x), x) is True
    assert is_function_smooth(sp.Abs(x), x) is False


def test_matrix_predicates_cover_definiteness_rank_and_parameter_strata() -> None:
    x, a = sp.symbols("x a", real=True)
    positive = sp.Matrix([[2, 1], [1, 2]])

    assert matrix_pd_on(positive) is True
    assert matrix_psd_on(positive) is True
    assert matrix_psd_on(sp.Matrix([[x]]), (x,), domain=x >= 0) is True
    assert matrix_rank_on(sp.Matrix([[1, x], [0, 0]]), (x,)) == 1
    ranks = matrix_rank_stratification(sp.Matrix([[a, 0], [0, 1]]), (a,))
    assert ranks.select({a: 0}) == 1
    assert ranks.select({a: 2}) == 2


def test_affine_relative_interior_and_strict_feasibility_share_the_hull() -> None:
    x, y = sp.symbols("x y", real=True)
    region = sp.And(sp.Eq(y, 0), x >= 0, x <= 1)
    relative = affine_relative_interior_formula(region, (x, y))
    result = strict_feasible(region, (x, y), return_result=True)

    assert relative.has(sp.Eq(y, 0))
    assert result.feasible
    assert result.witness is not None
    assert result.witness[y] == 0
    assert 0 < result.witness[x] < 1


def test_zero_dimensional_predicate_distinguishes_point_and_curve() -> None:
    x, y = sp.symbols("x y", real=True)

    assert is_zero_dimensional((x - 1, y + 2), (x, y)) is True
    assert is_zero_dimensional((x * y,), (x, y)) is False


def test_seeded_random_conveniences_are_reproducible_and_inside() -> None:
    ball = Ball((0, 0), 1)
    point = random_point(ball, seed=31)
    points = random_points(ball, count=3, seed=31)

    assert point == points[0]
    assert len(points) == 3
    assert all(ball.contains(candidate) for candidate in points)


def test_discretization_and_plotting_return_stable_public_objects() -> None:
    x, y = sp.symbols("x y", real=True)
    box = Box(((0, 1), (0, 1)))
    geometry = discretize_region_geometry(box, variables=(x, y))
    solution = SimpleNamespace(
        variables=(x, y),
        formula=(x >= 0) & (x <= 1) & (y >= 0) & (y <= 1),
        cells=(),
        samples=(),
    )

    assert geometry.dimension == 2
    assert len(geometry.polygons) == 1
    fig, (left, right) = plt.subplots(1, 2)
    try:
        assert plot_region_geometry(box, variables=(x, y), ax=left) is left
        assert plot_solution(solution, bounds=((0, 1), (0, 1)), ax=right) is right
        assert left.patches or left.collections or left.lines
        assert right.patches or right.collections or right.lines
    finally:
        plt.close(fig)
