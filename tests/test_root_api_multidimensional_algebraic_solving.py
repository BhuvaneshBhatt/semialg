"""Independent dimensions for multidimensional algebraic/solving APIs."""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib
import sympy as sp
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg import (
    Ball,
    Box,
    Interval,
    affine_relative_interior_formula,
    certified_radicalization,
    discretize_region_geometry,
    function_smoothness,
    irreducible_components,
    is_zero_dimensional,
    local_dimension_strata,
    matrix_definiteness,
    matrix_rank_stratification,
    minimal_prime_intersections,
    plot_region_geometry,
    plot_solution,
    random_point,
    random_points,
    reduced_component_singular_loci,
)

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt


def test_certified_radicalization_is_invariant_under_positive_multiplicity():
    x = sp.Symbol("x", real=True)
    assert (
        certified_radicalization(((x - 2) ** 2,), (x,)).equations
        == certified_radicalization(((x - 2) ** 5,), (x,)).equations
        == (x - 2,)
    )


def test_irreducible_components_ignore_nonzero_scalar_and_multiplicity():
    x, y = sp.symbols("x y", real=True)
    left = {c.equations for c in irreducible_components((7 * x**3 * y**2,), (x, y))}
    right = {c.equations for c in irreducible_components((x * y,), (x, y))}
    assert left == right == {(x,), (y,)}


def test_reduced_component_singular_loci_detect_intrinsic_cusp_not_crossing():
    x, y = sp.symbols("x y", real=True)
    loci = reduced_component_singular_loci((x * (y**2 - x**3),), (x, y))
    assert len(loci) == 2
    assert any(item.formula is not sp.false for item in loci)


def test_minimal_prime_intersections_respect_max_order():
    x, y, z = sp.symbols("x y z", real=True)
    # xyz=0 is the union of three coordinate planes.
    pairwise = minimal_prime_intersections((x * y * z,), (x, y, z), max_order=2)
    all_orders = minimal_prime_intersections((x * y * z,), (x, y, z))
    assert len(pairwise) == 3
    assert len(all_orders) == 4
    assert sorted(item.dimension for item in pairwise) == [1, 1, 1]
    assert all_orders[-1].dimension == 0


def test_local_dimension_strata_are_invariant_under_equation_scaling():
    x, y, z = sp.symbols("x y z", real=True)
    a = local_dimension_strata((x * y, x * z), (x, y, z))
    b = local_dimension_strata((3 * x * y, -5 * x * z), (x, y, z))
    assert [(s.dimension, sp.simplify_logic(s.formula)) for s in a] == [
        (s.dimension, sp.simplify_logic(s.formula)) for s in b
    ]


@given(st.integers(-4, 4), st.integers(-4, 4).filter(lambda n: n != 0))
@settings(max_examples=18, deadline=None)
def test_is_zero_dimensional_is_translation_and_scale_invariant(a, scale):
    x, y = sp.symbols("x y", real=True)
    assert is_zero_dimensional((scale * (x - a), y + a), (x, y)) is True
    assert is_zero_dimensional((scale * (x - a),), (x, y)) is False


def test_discretize_region_geometry_preserves_exact_box_vertices_under_translation():
    x, y = sp.symbols("x y", real=True)
    data = discretize_region_geometry(Box(((2, 5), (-3, 1))), variables=(x, y))
    assert set(data.polygons[0]) == {(2, -3), (5, -3), (5, 1), (2, 1)}


def test_random_point_is_exactly_first_seeded_random_points_sample():
    region = Interval(-3, 7)
    for seed in range(5):
        assert random_point(region, seed=seed) == random_points(region, count=4, seed=seed)[0]


@given(st.integers(0, 1000))
@settings(max_examples=12, deadline=None)
def test_random_points_seeded_samples_respect_ball_membership(seed):
    region = Ball((1, -2), 3)
    points = random_points(region, count=4, seed=seed)
    assert len(points) == 4
    assert all(region.contains(point) for point in points)


def test_plot_region_geometry_does_not_change_discretized_semantics():
    x, y = sp.symbols("x y", real=True)
    region = Box(((0, 2), (-1, 1)))
    before = discretize_region_geometry(region, variables=(x, y))
    fig, ax = plt.subplots()
    try:
        assert plot_region_geometry(region, variables=(x, y), ax=ax) is ax
    finally:
        plt.close(fig)
    after = discretize_region_geometry(region, variables=(x, y))
    assert before == after


def test_plot_solution_is_presentation_only_and_returns_supplied_axes():
    x = sp.Symbol("x", real=True)
    solution = SimpleNamespace(variables=(x,), formula=(x >= -1) & (x <= 1), cells=(), samples=())
    fig, ax = plt.subplots()
    try:
        assert plot_solution(solution, bounds=((-2, 2),), ax=ax) is ax
        assert ax.get_xlabel() == "x"
    finally:
        plt.close(fig)


def test_function_smoothness_is_representation_invariant_for_polynomials():
    x = sp.Symbol("x", real=True)
    a = function_smoothness((x - 1) ** 2, x)
    b = function_smoothness(x**2 - 2 * x + 1, x)
    assert (a.continuous, a.smooth, a.derivative_exceptions) == (
        b.continuous,
        b.smooth,
        b.derivative_exceptions,
    )
    assert a.smooth is True


def test_affine_relative_interior_formula_excludes_segment_endpoints_only():
    x, y = sp.symbols("x y", real=True)
    formula = affine_relative_interior_formula(sp.Eq(y, 0) & (x >= 0) & (x <= 2), (x, y))
    assert bool(formula.subs({x: 1, y: 0}))
    assert not bool(formula.subs({x: 0, y: 0}))
    assert not bool(formula.subs({x: 2, y: 0}))


def test_matrix_definiteness_is_congruence_invariant_for_constant_matrices():
    matrix = sp.Matrix([[2, 1], [1, 2]])
    transform = sp.Matrix([[1, 2], [0, 1]])
    assert matrix_definiteness(matrix, requested="positive_definite") is True
    assert (
        matrix_definiteness(transform.T * matrix * transform, requested="positive_definite") is True
    )


def test_matrix_rank_stratification_tracks_exact_rank_drop_boundary():
    a = sp.Symbol("a", real=True)
    result = matrix_rank_stratification(sp.diag(a**2, a - 1), (a,))
    assert result.select({a: 0}) == 1
    assert result.select({a: 1}) == 1
    assert result.select({a: 2}) == 2
