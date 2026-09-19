import pytest
import sympy as sp

from semialg import Circle, Simplex, Sphere


def test_sphere_through_planar_points_is_exact():
    sphere = Sphere.through((1, 0), (0, 1), (-1, 0))
    assert sphere == Sphere((0, 0), 1)
    assert all(sphere.contains(point) for point in ((1, 0), (0, 1), (-1, 0)))


def test_sphere_through_four_points_in_r3():
    sphere = Sphere.through((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
    assert sphere.center == (sp.Rational(1, 2),) * 3
    assert sphere.radius == sp.sqrt(3) / 2


def test_sphere_through_requires_unique_full_dimensional_data():
    with pytest.raises(ValueError, match="n \\+ 1"):
        Sphere.through((0, 0), (1, 0))
    with pytest.raises(ValueError, match="affinely independent"):
        Sphere.through((0, 0), (1, 0), (2, 0))


def test_circle_is_constructor_namespace_returning_sphere():
    circle = Circle((1, 2), 3)
    assert type(circle) is Sphere
    assert circle == Sphere((1, 2), 3)
    through = Circle.through((0, 0), (2, 0), (0, 2))
    assert type(through) is Sphere
    assert through.center == (1, 1)
    assert through.radius == sp.sqrt(2)


def test_circle_rejects_nonplanar_or_wrong_point_count():
    with pytest.raises(ValueError, match="two-dimensional center"):
        Circle((0, 0, 0), 1)
    with pytest.raises(ValueError, match="exactly three"):
        Circle.through((0, 0), (1, 0))
    with pytest.raises(ValueError, match="two-dimensional points"):
        Circle.through((0, 0, 0), (1, 0, 0), (0, 1, 0))


def test_circumscribed_triangle_matches_through_constructor():
    triangle = Simplex(((0, 0), (2, 0), (0, 2)))
    assert Sphere.circumscribed(triangle) == Circle.through(*triangle.vertices)


def test_inscribed_right_triangle_has_expected_center_and_radius():
    triangle = Simplex(((0, 0), (3, 0), (0, 4)))
    sphere = Sphere.inscribed(triangle)
    assert sphere == Sphere((1, 1), 1)


def test_inscribed_regular_tetrahedron_is_equidistant_from_faces():
    tetra = Simplex(((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)))
    sphere = Sphere.inscribed(tetra)
    assert sphere.center == (0, 0, 0)
    assert sp.simplify(sphere.radius - sp.sqrt(3) / 3) == 0


def test_simplex_sphere_constructors_require_full_dimension():
    embedded = Simplex(((0, 0, 0), (1, 0, 0), (0, 1, 0)))
    with pytest.raises(ValueError, match="full-dimensional"):
        Sphere.circumscribed(embedded)
    with pytest.raises(ValueError, match="full-dimensional"):
        Sphere.inscribed(embedded)
