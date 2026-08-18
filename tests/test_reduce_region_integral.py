import sympy as sp

from semialg import integrate_over_region, reduce_region_integral
from semialg.region_integrate import ReducedRegionIntegral


def test_reduce_region_integral_one_dimensional_intervals():
    x = sp.symbols("x", real=True)
    reduced = reduce_region_integral(x**2, x**2 <= 1, [x])
    assert isinstance(reduced, ReducedRegionIntegral)
    assert len(reduced.pieces) == 1
    assert reduced.pieces[0].limits == ((x, -1, 1),)
    assert sp.simplify(reduced.pieces[0].as_integral().doit() - sp.Rational(2, 3)) == 0


def test_reduce_region_integral_vertical_slice_triangle():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x >= 0, y >= 0, x + y <= 1)
    reduced = reduce_region_integral(x * y, condition, [x, y])
    assert isinstance(reduced, ReducedRegionIntegral)
    assert len(reduced.pieces) == 1
    piece = reduced.pieces[0]
    assert piece.limits == ((y, 0, 1 - x), (x, 0, 1))
    assert sp.simplify(piece.as_integral().doit() - sp.Rational(1, 24)) == 0


def test_reduce_region_integral_disk_as_vertical_slice():
    x, y = sp.symbols("x y", real=True)
    reduced = reduce_region_integral(1, x**2 + y**2 <= 1, [x, y])
    assert isinstance(reduced, ReducedRegionIntegral)
    assert len(reduced.pieces) == 1
    assert reduced.method == "radial_region_as_signed_vertical_slices"
    assert sp.simplify(reduced.unevaluated_sum().doit() - sp.pi) == 0


def test_reduce_region_integral_annulus_uses_signed_pieces():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x**2 + y**2 <= 4, x**2 + y**2 >= 1)
    reduced = reduce_region_integral(1, condition, [x, y])
    assert isinstance(reduced, ReducedRegionIntegral)
    assert len(reduced.pieces) == 2
    assert reduced.pieces[1].integrand == -1
    assert sp.simplify(reduced.unevaluated_sum().doit() - 3 * sp.pi) == 0


def test_integrate_over_region_delegates_to_reduction_layer():
    x, y = sp.symbols("x y", real=True)
    result = integrate_over_region(x**2 + y**2, x**2 + y**2 <= 1, [x, y], return_result=True)
    assert result.method == "radial_region_as_signed_vertical_slices"
    assert "pieces" in result.diagnostics
    assert sp.simplify(result.value - sp.pi / 2) == 0
