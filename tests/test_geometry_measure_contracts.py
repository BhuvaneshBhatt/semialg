import pytest
import sympy as sp

from semialg import (
    Ball,
    BoxRegion,
    IntervalRegion,
    Parallelepiped,
    Point,
    Polygon,
    Simplex,
    Sphere,
    SphericalShellRegion,
)
from semialg.standard_region_integrate import integrate_over_standard_region


@pytest.mark.parametrize(
    ("region", "variables"),
    [
        (Point((2, -1)), sp.symbols("x y", real=True)),
        (IntervalRegion(-2, 3), (sp.Symbol("x", real=True),)),
        (BoxRegion(((0, 2), (-1, 3))), sp.symbols("x y", real=True)),
        (Simplex(((0, 0), (2, 0), (0, 3))), sp.symbols("x y", real=True)),
        (Polygon(((0, 0), (2, 0), (2, 1), (0, 1))), sp.symbols("x y", real=True)),
        (Parallelepiped((0, 0), ((2, 0), (0, 3))), sp.symbols("x y", real=True)),
        (Ball((0, 0), 2), sp.symbols("x y", real=True)),
        (Sphere((0, 0), 2), sp.symbols("x y", real=True)),
        (SphericalShellRegion((0, 0), (1, 2)), sp.symbols("x y", real=True)),
    ],
    ids=[
        "point",
        "interval",
        "box",
        "simplex",
        "polygon",
        "parallelotope",
        "ball",
        "circle",
        "shell",
    ],
)
def test_measure_matches_integral(region, variables):
    direct = region.measure()
    integrated = integrate_over_standard_region(1, region, variables)
    assert sp.simplify(direct - integrated) == 0


@pytest.mark.parametrize(
    ("region", "variables"),
    [
        (IntervalRegion(-2, 4), (sp.Symbol("x", real=True),)),
        (BoxRegion(((0, 2), (-1, 3))), sp.symbols("x y", real=True)),
        (Simplex(((0, 0), (2, 0), (0, 3))), sp.symbols("x y", real=True)),
        (Simplex(((0, 0, 1), (2, 0, 1), (0, 3, 1))), sp.symbols("x y z", real=True)),
        (Polygon(((0, 0), (2, 0), (2, 1), (0, 1))), sp.symbols("x y", real=True)),
        (Sphere((1, -2), 3), sp.symbols("x y", real=True)),
    ],
    ids=["interval", "box", "simplex", "embedded-simplex", "polygon", "circle"],
)
def test_centroid_matches_moments(region, variables):
    measure = region.measure()
    expected = tuple(
        sp.simplify(integrate_over_standard_region(var, region, variables) / measure)
        for var in variables
    )
    assert region.centroid() == expected


@pytest.mark.parametrize(
    ("region", "intrinsic", "ambient"),
    [
        (Point((1, 2, 3)), 1, 0),
        (IntervalRegion(2, 2), 1, 0),
        (BoxRegion(((0, 0), (0, 2))), 2, 0),
        (BoxRegion(((1, 1), (2, 2))), 1, 0),
        (Simplex(((0, 0), (3, 4))), 5, 0),
        (Simplex(((0, 0, 1), (1, 0, 1), (0, 1, 1))), sp.Rational(1, 2), 0),
        (Sphere((0, 0), 2), 4 * sp.pi, 0),
        (Sphere((0, 0, 0), 2), 16 * sp.pi, 0),
        (Ball((0, 0, 0), 2), sp.Rational(32, 3) * sp.pi, sp.Rational(32, 3) * sp.pi),
    ],
    ids=[
        "point",
        "singleton",
        "segment-box",
        "point-box",
        "segment",
        "triangle",
        "circle",
        "sphere",
        "ball",
    ],
)
def test_dimension_measure_matrix(region, intrinsic, ambient):
    assert sp.simplify(region.measure() - intrinsic) == 0
    assert sp.simplify(region.measure(measure_dimension="ambient") - ambient) == 0
