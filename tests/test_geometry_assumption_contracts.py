import pytest
import sympy as sp

from semialg import Ball, Cylinder, Ellipsoid, SphericalShellRegion, Torus, Triangle


@pytest.mark.parametrize(
    ("region", "sufficient", "contradictory"),
    [
        (lambda r: Ball((0, 0), r), lambda r: r >= 0, lambda r: r < 0),
        (
            lambda r: SphericalShellRegion((0, 0), (0, r)),
            lambda r: r >= 0,
            lambda r: r < 0,
        ),
        (
            lambda r: Cylinder((0, 0, 0), (0, 0, 1), r),
            lambda r: r >= 0,
            lambda r: r < 0,
        ),
    ],
    ids=["ball", "shell", "cylinder"],
)
def test_radius_validity_pairs(region, sufficient, contradictory):
    r = sp.Symbol("r", real=True)
    geometry = region(r)
    assert geometry.is_valid() is None
    assert geometry.is_valid(sufficient(r)) is True
    assert geometry.is_valid(contradictory(r)) is False


def test_ellipsoid_validity_pair():
    a, b = sp.symbols("a b", real=True)
    region = Ellipsoid((0, 0), sp.diag(a, b))
    assert region.is_valid() is None
    assert region.is_valid(sp.And(a > 0, b > 0)) is True
    assert region.is_valid(a <= 0) is False


def test_triangle_validity_pair():
    a, b, c = sp.symbols("a b c", positive=True)
    region = Triangle.from_sides(a, b, c)
    sufficient = sp.And(a + b > c, a + c > b, b + c > a)
    contradictory = a + b <= c
    assert region.is_valid() is None
    assert region.is_valid(sufficient) is True
    assert region.is_valid(contradictory) is False


def test_torus_validity_pair():
    major, minor = sp.symbols("R r", real=True)
    region = Torus((0, 0, 0), major, minor)
    sufficient = sp.And(major > 0, minor > 0, major >= minor)
    contradictory = major < minor
    assert region.is_valid() is None
    assert region.is_valid(sufficient) is True
    assert region.is_valid(contradictory) is False
