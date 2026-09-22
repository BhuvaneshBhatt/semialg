import pytest
import sympy as sp

from semialg import (
    Ball,
    Box,
    Ellipsoid,
    EllipsoidBoundary,
    FinitePointSet,
    Interval,
    Parallelepiped,
    Polygon,
    Polytope,
    Simplex,
    Sphere,
    SphericalShell,
    TetrahedralComplex,
)
from semialg.standard_regions import Parallelogram


def _floats(point):
    return tuple(float(sp.N(value)) for value in point)


def _means(points):
    columns = tuple(zip(*(_floats(point) for point in points), strict=True))
    return tuple(sum(column) / len(column) for column in columns)


def _valid_sample(region, point):
    values = _floats(point)
    if isinstance(region, EllipsoidBoundary):
        center = tuple(float(sp.N(value)) for value in region.center)
        delta = sp.Matrix([value - origin for value, origin in zip(values, center, strict=True)])
        shape_inverse = sp.Matrix(region.shape_matrix).inv()
        level = float(sp.N((delta.T * shape_inverse * delta)[0]))
        return abs(level - 1.0) < 1e-10
    if isinstance(region, Sphere):
        center = tuple(float(sp.N(value)) for value in region.center)
        radius = float(sp.N(region.radius))
        squared = sum((value - origin) ** 2 for value, origin in zip(values, center, strict=True))
        return abs(squared - radius**2) < 1e-10
    if isinstance(region, Polytope):
        hrep = region.h_representation()
        vector = sp.Matrix(point)
        return all(
            float(sp.N((hrep.matrix.row(row) * vector)[0] - hrep.offsets[row, 0])) <= 1e-10
            for row in range(hrep.matrix.rows)
        )
    return bool(region.contains(point))


@pytest.mark.parametrize(
    ("region", "expected"),
    [
        (FinitePointSet(((0,), (2,))), (1.0,)),
        (Interval(0, 2), (1.0,)),
        (Box(((0, 2), (-1, 1))), (1.0, 0.0)),
        (Simplex(((0, 0), (1, 0), (0, 1))), (1 / 3, 1 / 3)),
        (Ball((1, -2), 2), (1.0, -2.0)),
        (Sphere((1, -2, 3), 2), (1.0, -2.0, 3.0)),
        (SphericalShell((1, -2), (1, 3)), (1.0, -2.0)),
        (Ellipsoid((1, -2), sp.diag(4, 9)), (1.0, -2.0)),
        (EllipsoidBoundary((1, -2), sp.diag(4, 9)), (1.0, -2.0)),
        (Polygon(((0, 0), (2, 0), (2, 1), (0, 1))), (1.0, 0.5)),
        (
            TetrahedralComplex((Simplex(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))),)),
            (0.25, 0.25, 0.25),
        ),
        (Parallelogram((0, 0), ((2, 0), (0, 4))), (1.0, 2.0)),
        (Parallelepiped((0, 0, 0), ((2, 0, 0), (0, 4, 0), (0, 0, 6))), (1.0, 2.0, 3.0)),
        (Polytope(((0, 0), (2, 0), (2, 2), (0, 2))), (1.0, 1.0)),
    ],
    ids=[
        "points",
        "interval",
        "box",
        "simplex",
        "ball",
        "sphere",
        "shell",
        "ellipsoid",
        "ellipsoid-boundary",
        "polygon",
        "polyhedron",
        "parallelogram",
        "parallelepiped",
        "polytope",
    ],
)
def test_uniform_sampling_mean(region, expected):
    points = region.random_points(1000, seed=2718)
    assert all(_valid_sample(region, point) for point in points[:6])
    for observed, target in zip(_means(points), expected, strict=True):
        assert abs(observed - target) < 0.13


def test_seed_reproducibility():
    region = Ball((0, 0), 1)
    assert region.random_points(12, seed=91) == region.random_points(12, seed=91)


def test_ball_radial_moment():
    points = Ball((0, 0), 1).random_points(4000, seed=7)
    mean_radius_squared = sum(sum(value * value for value in _floats(p)) for p in points) / len(
        points
    )
    assert abs(mean_radius_squared - 0.5) < 0.025


def test_shell_radial_moment():
    points = SphericalShell((0, 0), (1, 3)).random_points(4000, seed=17)
    mean_radius_squared = sum(sum(value * value for value in _floats(p)) for p in points) / len(
        points
    )
    assert abs(mean_radius_squared - 5.0) < 0.12


def test_sphere_second_moment():
    points = Sphere((0, 0, 0), 2).random_points(4000, seed=31)
    mean_x_squared = sum(_floats(point)[0] ** 2 for point in points) / len(points)
    assert abs(mean_x_squared - 4 / 3) < 0.08


def test_ellipse_surface_moment():
    ellipse = EllipsoidBoundary((0, 0), sp.diag(4, 9))
    points = ellipse.random_points(5000, seed=41)
    mean_x_squared = sum(_floats(point)[0] ** 2 for point in points) / len(points)
    # Arc-length weighting gives about 2.19702. Uniform angular sampling would
    # give 2, so this contract detects the common affine-sphere bias.
    assert abs(mean_x_squared - 2.19702) < 0.08


def test_degenerate_parallelotope_declines():
    region = Parallelogram((0, 0), ((1, 0), (2, 0)))
    with pytest.raises(NotImplementedError, match="independent spanning vectors"):
        region.random_point(seed=1)
