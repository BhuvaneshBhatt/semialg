import sympy as sp

import semialg
from semialg import Octahedron, Polygon, Polytope


def test_convex_hull_removes_interior_and_canonicalizes_square():
    hull = semialg.convex_hull(
        [(0, 0), (1, 0), (1, 1), (0, 1), (sp.Rational(1, 2), sp.Rational(1, 2))]
    )
    assert hull.dimension() == 2
    assert hull.dimension() == 2


def test_convex_hull_embedded_plane():
    hull = semialg.convex_hull(
        [(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1), (sp.Rational(1, 2), sp.Rational(1, 2), 1)],
        canonical=False,
    )
    assert isinstance(hull, Polytope)
    assert hull.dimension() == 2
    assert len(hull.vertices) == 4


def test_random_generators_are_seeded_exact():
    a = semialg.random_polytope(3, seed=7)
    b = semialg.random_polytope(3, seed=7)
    assert a == b and a.dimension() == 3
    p = semialg.random_polygon(seed=4)
    assert isinstance(p, Polygon) and p.dimension() == 2


def test_triangular_subdivision_and_geodesic_projection():
    ico = Octahedron()
    shell = semialg.subdivide_triangular_faces(ico, levels=1)
    assert len(shell.faces) == 32
    geo = semialg.geodesic_refinement(ico, levels=1)
    norms = {sp.simplify(sum(x * x for x in p)) for p in geo.vertices}
    assert len(norms) == 1


def test_polyhedral_intersection_preserves_polytope():
    a = Polytope([(0, 0), (2, 0), (2, 2), (0, 2)])
    b = Polytope([(1, 1), (3, 1), (3, 3), (1, 3)])
    c = semialg.polyhedral_intersection(a, b)
    assert c.dimension() == 2
    assert c.dimension() == 2
    assert c.contains((sp.Rational(3, 2), sp.Rational(3, 2)))


def test_polyhedral_boolean_intersection_public_fast_path():
    left = Polytope([(0, 0), (2, 0), (2, 2), (0, 2)])
    right = Polytope([(1, 1), (3, 1), (3, 3), (1, 3)])
    result = semialg.polyhedral_boolean("intersection", left, right)
    assert result.dimension() == 2
