"""Public API coverage for geometry generation and exact convex hull."""

import sympy as sp

import semialg
from semialg import Octahedron, Polytope


def test_exact_hull_and_seeded_generation_public_calls():
    hull = semialg.convex_hull([(0, 0), (2, 0), (0, 2), (1, 1)])
    assert hull.dimension() == 2
    assert semialg.random_polygon(seed=11).dimension() == 2
    assert semialg.random_polytope(2, seed=11).dimension() == 2


def test_refinement_and_polyhedral_boolean_public_calls():
    shell = semialg.subdivide_triangular_faces(Octahedron(), levels=1)
    assert len(shell.faces) == 32
    assert semialg.geodesic_refinement(Octahedron(), levels=1).dimension() == 2
    left = Polytope([(0, 0), (2, 0), (2, 2), (0, 2)])
    right = Polytope([(1, 1), (3, 1), (3, 3), (1, 3)])
    intersection = semialg.polyhedral_intersection(left, right)
    assert intersection.contains((sp.Rational(3, 2), sp.Rational(3, 2)))
    assert semialg.polyhedral_boolean("intersection", left, right).dimension() == 2
