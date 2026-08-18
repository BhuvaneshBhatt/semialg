import pytest
import sympy as sp

from semialg import integrate_over_region, semialgebraic_measure

pytestmark = pytest.mark.slow


def test_zero_dimensional_measure_counts_univariate_real_roots():
    x = sp.symbols("x", real=True)

    assert semialgebraic_measure(sp.Eq(x**2 - 1, 0), [x], measure_dimension=0) == 2
    assert integrate_over_region(x**2, sp.Eq(x**2 - 1, 0), [x], measure_dimension=0) == 2


def test_zero_dimensional_measure_counts_plane_point_system():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(sp.Eq(x, 0), sp.Eq(y, 2))

    assert semialgebraic_measure(condition, [x, y], measure_dimension=0) == 1
    assert integrate_over_region(x + y, condition, [x, y], measure_dimension=0) == 2


def test_intrinsic_circle_length_and_moment():
    x, y = sp.symbols("x y", real=True)
    circle = sp.Eq(x**2 + y**2, 1)

    assert (
        sp.simplify(
            semialgebraic_measure(circle, [x, y], measure_dimension="intrinsic") - 2 * sp.pi
        )
        == 0
    )
    assert (
        sp.simplify(integrate_over_region(x**2, circle, [x, y], measure_dimension=1) - sp.pi) == 0
    )


def test_graph_curve_intrinsic_length():
    x, y = sp.symbols("x y", real=True)
    segment = sp.And(sp.Eq(y, x), x >= 0, x <= 1)

    assert (
        sp.simplify(
            semialgebraic_measure(segment, [x, y], measure_dimension="intrinsic") - sp.sqrt(2)
        )
        == 0
    )
    assert (
        sp.simplify(integrate_over_region(x, segment, [x, y], measure_dimension=1) - sp.sqrt(2) / 2)
        == 0
    )


def test_ambient_measure_remains_default_for_lower_dimensional_sets():
    x, y = sp.symbols("x y", real=True)

    assert semialgebraic_measure(sp.Eq(x**2 + y**2, 1), [x, y]) == 0
