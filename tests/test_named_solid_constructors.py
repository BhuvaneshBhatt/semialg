import itertools

import pytest
import sympy as sp

from semialg import (
    Cube,
    Dodecahedron,
    Hexahedron,
    Icosahedron,
    Octahedron,
    Parallelepiped,
    Polytope,
    Prism,
    Pyramid,
    Simplex,
    Tetrahedron,
)


def _minimum_pair_distance(vertices):
    distances = [
        sp.simplify(sp.sqrt(sum((a - b) ** 2 for a, b in zip(p, q, strict=True))))
        for p, q in itertools.combinations(vertices, 2)
    ]
    return min(distances, key=lambda value: float(sp.N(value)))


def test_cube_returns_canonical_parallelepiped():
    cube = Cube((1, 2, 3), 2)
    assert isinstance(cube, Parallelepiped)
    assert cube.origin == (0, 1, 2)
    assert cube.vectors == ((2, 0, 0), (0, 2, 0), (0, 0, 2))
    assert cube.dimension() == 3


def test_regular_tetrahedron_has_requested_edge():
    tet = Tetrahedron(edge=sp.sqrt(2))
    assert isinstance(tet, Simplex)
    assert len(tet.vertices) == 4
    assert sp.simplify(_minimum_pair_distance(tet.vertices) - sp.sqrt(2)) == 0


def test_explicit_tetrahedron_preserves_vertices():
    vertices = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
    assert Tetrahedron(vertices).vertices == vertices


def _assert_regular_polytope(solid, count):
    assert isinstance(solid, Polytope)
    assert len(solid.vertices) == count
    assert solid.dimension() == 3
    assert sp.simplify(_minimum_pair_distance(solid.vertices) - 3) == 0


def test_regular_platonic_polytope_edges():
    _assert_regular_polytope(Octahedron(edge=3), 6)
    _assert_regular_polytope(Icosahedron(edge=3), 12)
    _assert_regular_polytope(Dodecahedron(edge=3), 20)


def test_prism_extrudes_base_vertices():
    prism = Prism(((0, 0, 0), (1, 0, 0), (0, 1, 0)), (0, 0, 2))
    assert isinstance(prism, Polytope)
    assert len(prism.vertices) == 6
    assert prism.dimension() == 3
    assert (0, 0, 2) in prism.vertices


def test_pyramid_requires_apex_outside_base_hull():
    pyramid = Pyramid(((0, 0, 0), (1, 0, 0), (0, 1, 0)), (0, 0, 2))
    assert isinstance(pyramid, Polytope)
    assert len(pyramid.vertices) == 4
    assert pyramid.dimension() == 3
    with pytest.raises(ValueError):
        Pyramid(((0, 0, 0), (1, 0, 0), (0, 1, 0)), (1, 1, 0))


def test_hexahedron_requires_eight_distinct_full_dimensional_vertices():
    cube_vertices = tuple(itertools.product((0, 1), repeat=3))
    hexa = Hexahedron(cube_vertices)
    assert isinstance(hexa, Polytope)
    assert len(hexa.vertices) == 8
    with pytest.raises(ValueError):
        Hexahedron(cube_vertices[:-1])


def test_named_solid_positive_lengths_are_validated():
    with pytest.raises(ValueError):
        Cube(side=0)
    with pytest.raises(ValueError):
        Icosahedron(edge=-1)
