import math

import sympy as sp

from semialg import (
    Ball,
    Box,
    Ellipsoid,
    Simplex,
    Sphere,
    SphericalShell,
    as_semialgebraic_region,
    random_point,
    random_points,
    region_measure,
)


def _float_tuple(point):
    return tuple(float(sp.N(value)) for value in point)


def test_geometry_representative_sample_is_tuple_and_valid():
    ball = Ball((0, 0), 2)
    point = ball.sample_point()
    assert isinstance(point, tuple)
    assert ball.contains(point)


def test_ball_random_point_is_seeded_and_inside():
    ball = Ball((1, -2), 3)
    left = ball.random_point(seed=17)
    right = ball.random_point(seed=17)
    assert left == right
    x, y = _float_tuple(left)
    assert (x - 1) ** 2 + (y + 2) ** 2 <= 9 + 1e-12


def test_sphere_random_point_uses_surface_measure():
    sphere = Sphere((0, 0, 0), 2)
    point = sphere.random_point(seed=11)
    norm = math.sqrt(sum(value * value for value in _float_tuple(point)))
    assert abs(norm - 2.0) < 1e-12


def test_shell_random_point_respects_shell():
    shell = SphericalShell((0, 0), (1, 3))
    for point in shell.random_points(8, seed=4):
        radius = math.hypot(*_float_tuple(point))
        assert 1 <= radius <= 3


def test_simplex_uniform_sampler_stays_in_simplex():
    simplex = Simplex(((0, 0), (1, 0), (0, 1)))
    for point in simplex.random_points(12, seed=3):
        x, y = _float_tuple(point)
        assert x >= 0 and y >= 0 and x + y <= 1 + 1e-12


def test_ellipsoid_random_point_and_measure():
    ellipsoid = Ellipsoid((2, -1), ((4, 0), (0, 9)))
    assert ellipsoid.measure() == 6 * sp.pi
    assert ellipsoid.centroid() == (2, -1)
    for point in ellipsoid.random_points(5, seed=9):
        x, y = _float_tuple(point)
        assert ((x - 2) ** 2) / 4 + ((y + 1) ** 2) / 9 <= 1 + 1e-12


def test_intrinsic_vs_ambient_measure_for_sphere():
    sphere = Sphere((0, 0), 2)
    assert sphere.measure() == 4 * sp.pi
    assert sphere.measure(measure_dimension="ambient") == 0


def test_box_measure_and_centroid():
    box = Box(((0, 2), (-1, 3)))
    assert region_measure(box) == 8
    assert box.centroid() == (1, 1)


def test_formula_region_random_sampling_with_explicit_bounds():
    x, y = sp.symbols("x y", real=True)
    region = as_semialgebraic_region(x**2 + y**2 <= 1, (x, y))
    point = region.random_point(seed=5, bounds=((-1, 1), (-1, 1)))
    px, py = _float_tuple(point)
    assert px * px + py * py <= 1 + 1e-12


def test_public_random_points_function_for_canonical_region():
    points = random_points(Ball((0, 0), 1), count=3, seed=2)
    assert len(points) == 3
    assert random_point(Ball((0, 0), 1), seed=2) == points[0]


def test_region_measure_default_depends_on_input_abstraction():
    x, y = sp.symbols("x y", real=True)
    segment = sp.And(sp.Eq(y, 0), x >= 0, x <= 1)

    assert region_measure(segment, (x, y)) == 0
    assert region_measure(segment, (x, y), measure_dimension="intrinsic") == 1

    circle = Sphere((0, 0), 2)
    assert region_measure(circle) == 4 * sp.pi
    assert region_measure(circle, measure_dimension="ambient") == 0
