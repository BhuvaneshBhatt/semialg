from __future__ import annotations

import sympy as sp

import semialg.region_integrate  # noqa: F401 - initialize integration module graph
from semialg._region_integrate_geometry import _vertical_slice_data
from semialg._region_integrate_intrinsic import _circle_radius_squared


def test_vertical_slice_declines_undecidable_parameter_equality() -> None:
    x, y, a = sp.symbols("x y a", real=True)
    condition = sp.And(sp.Eq(a, 0), y >= 0, y <= 1, x >= 0, x <= 1)
    assert _vertical_slice_data(condition, x, y, {}) is None


def test_circle_recognizer_requires_certified_matching_quadratic_coefficients() -> None:
    x, y, a = sp.symbols("x y a", real=True)
    condition = sp.Eq(x**2 + a * y**2 - 1, 0)
    assert _circle_radius_squared(condition, x, y) is None


def test_vertical_slice_preserves_base_variable_equality() -> None:
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(sp.Eq(x, 0), y >= 0, y <= 1)
    data = _vertical_slice_data(condition, x, y, {})
    assert data is not None
    lower, upper, intervals = data
    assert (lower, upper) == (0, 1)
    assert intervals == ()
