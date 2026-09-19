import itertools

import pytest
import sympy as sp

from semialg import Hexahedron, HRepresentation, Polytope


def test_cube_v_to_h_and_incidence_are_exact():
    cube = Polytope(tuple(itertools.product((0, 1), repeat=3)))
    hrep = cube.h_representation()
    assert isinstance(hrep, HRepresentation)
    assert hrep.matrix.shape == (6, 3)
    incidence = cube.incidence()
    assert incidence.facet_count == 6
    assert sorted(len(f.vertex_indices) for f in incidence.facets) == [4] * 6
    assert all(len(adjacent) == 3 for adjacent in incidence.vertex_facets)


def test_h_to_v_round_trip_for_box():
    A = ((1, 0), (-1, 0), (0, 1), (0, -1))
    b = (2, 0, 3, 0)
    poly = Polytope.from_halfspaces(A, b)
    assert set(poly.vertices) == {(0, 0), (2, 0), (0, 3), (2, 3)}
    assert poly.validate_topology(vertex_count=4, facet_count=4, facet_sizes=(2,) * 4)


def test_unbounded_halfspaces_are_rejected():
    with pytest.raises(ValueError, match="bounded"):
        Polytope.from_halfspaces(((1, 0), (-1, 0)), (1, 0))


def test_exact_algebraic_h_representation_round_trip():
    r = sp.sqrt(2)
    poly = Polytope.from_halfspaces(((1, 0), (-1, 0), (0, 1), (0, -1)), (r, 0, r, 0))
    assert (r, r) in poly.vertices
    assert set(poly.h_representation().vertices()) == set(poly.vertices)


def test_hexahedron_rejects_eight_point_non_hexahedral_hull():
    cube = list(itertools.product((0, 1), repeat=3))
    cube[-1] = (2, 2, 2)
    with pytest.raises(ValueError, match="six quadrilateral facets"):
        Hexahedron(cube)
