from __future__ import annotations

import pytest
import sympy as sp

from semialg import region_centroid, region_covariance, region_moment

pytestmark = pytest.mark.slow


def test_raw_moments_on_interval() -> None:
    x = sp.symbols("x", real=True)

    assert sp.simplify(region_moment(x**2 <= 1, [x], powers=[2]) - sp.Rational(2, 3)) == 0
    assert sp.simplify(region_moment(x**2 <= 1, [x]) - 2) == 0


def test_centroid_and_covariance_on_symmetric_interval() -> None:
    x = sp.symbols("x", real=True)

    centroid = region_centroid(x**2 <= 1, [x])
    assert centroid[x] == 0

    covariance = region_covariance(x**2 <= 1, [x])
    assert covariance == sp.Matrix([[sp.Rational(1, 3)]])


def test_centroid_and_covariance_on_unit_square() -> None:
    x, y = sp.symbols("x y", real=True)
    square = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)

    centroid = region_centroid(square, [x, y])
    assert centroid == {x: sp.Rational(1, 2), y: sp.Rational(1, 2)}

    covariance = region_covariance(square, [x, y])
    assert covariance == sp.Matrix([[sp.Rational(1, 12), 0], [0, sp.Rational(1, 12)]])


def test_centroid_on_standard_simplex() -> None:
    x, y = sp.symbols("x y", real=True)
    triangle = sp.And(x >= 0, y >= 0, x + y <= 1)

    centroid = region_centroid(triangle, [x, y])
    assert centroid == {x: sp.Rational(1, 3), y: sp.Rational(1, 3)}

    covariance = region_covariance(triangle, [x, y])
    assert covariance == sp.Matrix(
        [[sp.Rational(1, 18), -sp.Rational(1, 36)], [-sp.Rational(1, 36), sp.Rational(1, 18)]]
    )


def test_disk_moments_centroid_covariance() -> None:
    x, y = sp.symbols("x y", real=True)
    disk = x**2 + y**2 <= 1

    assert sp.simplify(region_moment(disk, [x, y], powers=[2, 0]) - sp.pi / 4) == 0
    assert region_centroid(disk, [x, y]) == {x: 0, y: 0}
    assert region_covariance(disk, [x, y]) == sp.Matrix(
        [[sp.Rational(1, 4), 0], [0, sp.Rational(1, 4)]]
    )


def test_intrinsic_circle_centroid_and_covariance() -> None:
    x, y = sp.symbols("x y", real=True)
    circle = sp.Eq(x**2 + y**2, 1)

    centroid = region_centroid(circle, [x, y], measure_dimension="intrinsic")
    assert centroid == {x: 0, y: 0}

    covariance = region_covariance(circle, [x, y], measure_dimension=1)
    assert covariance == sp.Matrix([[sp.Rational(1, 2), 0], [0, sp.Rational(1, 2)]])
