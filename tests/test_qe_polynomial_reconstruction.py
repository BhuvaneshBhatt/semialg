"""Polynomial reconstruction regressions for CAD-derived formulas."""

from __future__ import annotations

import sympy as sp

from semialg import (
    is_equal,
    semialgebraic_image,
    semialgebraic_projection,
)


def test_linear_image_reconstructs_as_polynomial_disk_and_reenters_cad():
    x, y, u, v = sp.symbols("x y u v", real=True)

    image = semialgebraic_image(
        [x + y, x - y],
        x**2 + y**2 <= 1,
        [x, y],
        image_variables=[u, v],
    )
    expected = u**2 + v**2 <= 2

    assert image == (u**2 + v**2 - 2 <= 0)
    assert is_equal(image, expected, [u, v])


def test_projection_recovers_single_polynomial_sign_relation():
    x, y = sp.symbols("x y", real=True)

    projected = semialgebraic_projection(
        (y >= x**2) & (y <= 1),
        eliminate=[y],
        variables=[x, y],
    )

    assert is_equal(projected, x**2 <= 1, [x])


def test_polynomial_sign_reconstruction_preserves_asymmetric_ray():
    x, y = sp.symbols("x y", real=True)

    projected = semialgebraic_projection(
        sp.Eq(y, x) & (x >= 1),
        eliminate=[y],
        variables=[x, y],
    )

    assert is_equal(projected, x >= 1, [x])
    assert sp.simplify(projected.subs(x, -2)) is sp.false
