import pytest
import sympy as sp

from semialg import RegionIntegralResult, integrate_over_region

pytestmark = pytest.mark.slow


def test_integrate_over_one_dimensional_semialgebraic_set():
    x = sp.symbols("x", real=True)
    assert sp.simplify(integrate_over_region(x**2, x**2 <= 1, [x]) - sp.Rational(2, 3)) == 0
    assert sp.simplify(integrate_over_region(x, (x > 0) & (x < 1), [x]) - sp.Rational(1, 2)) == 0


def test_integrate_over_one_dimensional_boolean_union_with_bounds():
    x = sp.symbols("x", real=True)
    value = integrate_over_region(1 + x, (x < -1) | (x > 1), [x], bounds=[(x, -2, 2)])
    assert sp.simplify(value - 2) == 0


def test_integrate_over_unit_disk_polynomial_moments():
    x, y = sp.symbols("x y", real=True)
    assert sp.simplify(integrate_over_region(1, x**2 + y**2 <= 1, [x, y]) - sp.pi) == 0
    assert sp.simplify(integrate_over_region(x, x**2 + y**2 <= 1, [x, y])) == 0
    assert (
        sp.simplify(integrate_over_region(x**2 + y**2, x**2 + y**2 <= 1, [x, y]) - sp.pi / 2) == 0
    )


def test_integrate_over_annulus_polynomial_moments():
    x, y = sp.symbols("x y", real=True)
    condition = (x**2 + y**2 <= 4) & (x**2 + y**2 >= 1)
    assert sp.simplify(integrate_over_region(1, condition, [x, y]) - 3 * sp.pi) == 0
    assert (
        sp.simplify(
            integrate_over_region(x**2 + y**2, condition, [x, y]) - sp.Rational(15, 2) * sp.pi
        )
        == 0
    )


def test_integrate_over_triangle_vertical_slice():
    x, y = sp.symbols("x y", real=True)
    condition = (x >= 0) & (y >= 0) & (x + y <= 1)
    assert sp.simplify(integrate_over_region(1, condition, [x, y]) - sp.Rational(1, 2)) == 0
    assert sp.simplify(integrate_over_region(x * y, condition, [x, y]) - sp.Rational(1, 24)) == 0


def test_integrate_over_parabola_cap_vertical_slice():
    x, y = sp.symbols("x y", real=True)
    condition = (y >= x**2) & (y <= 1)
    assert sp.simplify(integrate_over_region(1, condition, [x, y]) - sp.Rational(4, 3)) == 0
    assert sp.simplify(integrate_over_region(y, condition, [x, y]) - sp.Rational(4, 5)) == 0


def test_integrate_over_region_result_object():
    x = sp.symbols("x", real=True)
    result = integrate_over_region(x, (x >= 0) & (x <= 2), [x], return_result=True)
    assert isinstance(result, RegionIntegralResult)
    assert result.method == "one_dimensional_cell_integration"
    assert sp.simplify(result.value - 2) == 0
