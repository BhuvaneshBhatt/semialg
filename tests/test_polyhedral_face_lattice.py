import itertools

import sympy as sp

from semialg import HRepresentation, Polytope
from semialg.polyhedral import verify_h_redundancy_certificate


def _cube():
    return Polytope(tuple(itertools.product((0, 1), repeat=3)))


def test_cube_face_lattice_edges_ridges_and_f_vector():
    cube = _cube()
    assert cube.f_vector() == (8, 12, 6)
    assert len(cube.edges()) == 12
    assert len(cube.ridges()) == 12
    assert cube.euler_characteristic() == 1
    assert cube.face_lattice().boundary_euler_characteristic == 2
    assert cube.face_lattice().extended_f_vector == (1, 8, 12, 6, 1)


def test_cube_vertex_and_facet_adjacency():
    adjacency = _cube().adjacency()
    assert all(len(neighbors) == 3 for neighbors in adjacency.vertex_neighbors)
    assert all(len(neighbors) == 4 for neighbors in adjacency.facet_neighbors)


def test_arbitrary_k_faces_for_four_simplex():
    simplex = Polytope(((0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)))
    assert simplex.f_vector() == (5, 10, 10, 5)
    assert len(simplex.faces(2)) == 10
    assert len(simplex.ridges()) == 10


def test_redundancy_certificates_are_replayable():
    hrep = HRepresentation(
        ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 0)),
        (1, 0, 1, 0, 2),
    )
    certificates = hrep.redundancy_certificates()
    assert [c.redundant for c in certificates] == [False, False, False, False, True]
    assert all(verify_h_redundancy_certificate(hrep, cert) for cert in certificates)
    reduced = hrep.irredundant()
    assert reduced.matrix.rows == 4
    assert set(reduced.vertices()) == {(0, 0), (0, 1), (1, 0), (1, 1)}


def test_scaled_duplicate_constraint_is_redundant():
    hrep = HRepresentation(((1,), (-1,), (2,)), (1, 0, 2))
    certs = hrep.redundancy_certificates()
    assert sum(cert.redundant for cert in certs) == 2
    assert hrep.irredundant().matrix.rows == 2


def test_combinatorial_equivalence_ignores_affine_geometry():
    cube = _cube()
    skew = Polytope(
        tuple((2 * x + y, 3 * y + z, z) for x, y, z in itertools.product((0, 1), repeat=3))
    )
    assert cube.is_combinatorially_equivalent(skew)


def test_non_equivalent_polytopes_are_distinguished():
    cube = _cube()
    octahedron = Polytope(((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)))
    assert not cube.is_combinatorially_equivalent(octahedron)


def test_algebraic_coordinates_preserve_face_lattice():
    r = sp.sqrt(2)
    box = Polytope(tuple(itertools.product((0, r), repeat=3)))
    assert box.f_vector() == (8, 12, 6)
