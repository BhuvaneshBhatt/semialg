from semialg import (
    Box,
    Polytope,
    Simplex,
    canonicalize_polygon,
    canonicalize_polyhedron,
    canonicalize_region,
    decompose_polytope,
    tetrahedralize_cell,
    tetrahedralize_cells,
    triangulate_polytope,
)


def cube():
    return [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]


def test_polygon_triangle_canonicalizes_to_simplex():
    assert isinstance(canonicalize_polygon([(0, 0), (1, 0), (0, 1)]), Simplex)


def test_cube_canonicalizes_to_box():
    assert isinstance(canonicalize_polyhedron(Polytope(cube())), Box)


def test_pulling_cube_is_exact_tetrahedralization():
    d = triangulate_polytope(Polytope(cube()))
    assert d.dimension == 3
    assert len(d.simplices) == 6
    assert all(len(s) == 4 for s in d.simplices)


def test_placing_and_barycentric_available():
    p = Polytope([(0, 0), (1, 0), (1, 1), (0, 1)])
    assert len(triangulate_polytope(p, strategy="placing").simplices) == 2
    b = triangulate_polytope(p, strategy="barycentric")
    assert len(b.simplices) == 8
    assert b.introduced_vertex_indices


def test_decomposition_targets():
    assert decompose_polytope(Polytope(cube()), target="tetrahedra").dimension == 3


def test_mixed_cell_cube_count_and_global_ids():
    r = tetrahedralize_cell(cube(), vertex_ids=(10, 11, 12, 13, 14, 15, 16, 17))
    assert len(r.tetrahedra) == 6
    assert set(sum((tuple(t) for t in r.tetrahedra), ())) <= set(range(10, 18))


def test_adjacent_cubes_share_same_face_triangulation():
    pts = [
        (0, 0, 0),
        (1, 0, 0),
        (0, 1, 0),
        (1, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (0, 1, 1),
        (1, 1, 1),
        (2, 0, 0),
        (2, 1, 0),
        (2, 0, 1),
        (2, 1, 1),
    ]
    cells = [(0, 1, 2, 3, 4, 5, 6, 7), (1, 8, 3, 9, 5, 10, 7, 11)]
    r = tetrahedralize_cells(pts, cells)
    assert len(r.tetrahedra) == 12
    shared = {1, 3, 5, 7}
    tris = []
    for a, b in r.source_cell_ranges:
        boundary = []
        for t in r.tetrahedra[a:b]:
            for face in __import__("itertools").combinations(t, 3):
                if set(face) <= shared:
                    boundary.append(frozenset(face))
        tris.append(set(boundary))
    assert tris[0] == tris[1] and len(tris[0]) == 2


def test_simplex_pulling_is_one_simplex():
    p = Polytope([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
    assert len(triangulate_polytope(p).simplices) == 1


def test_triangular_prism_pulling_has_three_tetrahedra():
    p = Polytope([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)])
    assert len(triangulate_polytope(p).simplices) == 3


def test_square_pyramid_pulling_has_two_tetrahedra():
    p = Polytope([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1)])
    assert len(triangulate_polytope(p).simplices) == 2


def test_canonicalize_region_dispatch():
    assert isinstance(canonicalize_region([(0, 0), (1, 0), (0, 1)]), list)
