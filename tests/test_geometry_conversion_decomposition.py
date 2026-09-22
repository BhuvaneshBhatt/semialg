import sympy as sp

from semialg import (
    canonicalize_polygon,
    canonicalize_polyhedron,
    convert_region,
    polygonal_region_from_paths,
)
from semialg.boundary_topology import PolygonalSet, PolyhedralComponent, PolyhedralShell, Polyhedron
from semialg.mixed_tetrahedralization import tetrahedralize_cell, tetrahedralize_cells
from semialg.polytope_decomposition import decompose_polytope, triangulate_polytope
from semialg.standard_regions import Box, Parallelogram, Polygon, Polytope, Simplex


def cube(x0=0, x1=1):
    return [
        (x0, 0, 0),
        (x1, 0, 0),
        (x0, 1, 0),
        (x1, 1, 0),
        (x0, 0, 1),
        (x1, 0, 1),
        (x0, 1, 1),
        (x1, 1, 1),
    ]


def cube_shell():
    v = cube()
    f = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    return PolyhedralShell(v, f)


def test_p2_polygon_recognition():
    assert isinstance(canonicalize_polygon([(0, 0), (2, 0), (2, 1), (0, 1)]), Box)
    assert isinstance(canonicalize_polygon([(0, 0), (2, 0), (3, 1), (1, 1)]), Parallelogram)
    assert isinstance(canonicalize_polygon([(0, 0), (1, 0), (0, 1)]), Simplex)


def test_p2_nonconvex_shell_preserved():
    # Two tetrahedra represented as a nonconvex tetrahedral complex remain nonconvex representation.
    from semialg.standard_regions import TetrahedralComplex

    tc = TetrahedralComplex(
        [
            [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
            [(0, 0, 0), (-1, 0, 0), (0, 1, 0), (0, 0, 1)],
        ]
    )
    assert canonicalize_polyhedron(tc) is tc


def test_p3_crossing_and_fill_rules():
    bow = [(0, 0), (2, 2), (0, 2), (2, 0)]
    r = polygonal_region_from_paths([bow], fill_rule="nonzero")
    assert isinstance(r, PolygonalSet) and len(r.components) == 2
    doubled = [(0, 0), (2, 0), (2, 2), (0, 2)]
    r2 = polygonal_region_from_paths([doubled, doubled], fill_rule="winding_at_least_two")
    assert len(r2.components) == 1
    assert (
        len(polygonal_region_from_paths([doubled, doubled], fill_rule="even_odd").components) == 0
    )


def test_p4_conversion_contract():
    p = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    assert convert_region(p, "canonical") is not None
    assert convert_region(p, "formula") is not None
    assert convert_region(p, "semialgebraic").formula is not None
    assert convert_region(p, "boundary").components
    rr = convert_region(p, "formula", return_result=True)
    assert rr.exact and rr.certified and not rr.lossy


def test_p5_strategies_and_targets():
    p = Polytope(cube())
    for strategy in ("pulling", "placing", "barycentric"):
        d = triangulate_polytope(p, strategy=strategy)
        assert d.certified and d.conforming and all(len(s) == 4 for s in d.simplices)
    assert decompose_polytope(p, target="tetrahedra").dimension == 3


def test_p5_polyhedron_formula_convex_shell():
    x, y, z = sp.symbols("x y z", real=True)
    poly = Polyhedron((PolyhedralComponent(cube_shell()),))
    f = convert_region(poly, "formula", variables=(x, y, z))
    assert f is not sp.false


def test_p6_cell_counts_and_shared_face_conformity():
    tet = tetrahedralize_cell([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
    pyr = tetrahedralize_cell([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1)])
    prism = tetrahedralize_cell([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)])
    hexa = tetrahedralize_cell(cube())
    assert (
        len(tet.tetrahedra) == 1
        and len(pyr.tetrahedra) == 2
        and len(prism.tetrahedra) == 3
        and len(hexa.tetrahedra) in (5, 6)
    )
    verts = cube() + [(2, 0, 0), (2, 1, 0), (2, 0, 1), (2, 1, 1)]
    m = tetrahedralize_cells(verts, [(0, 1, 2, 3, 4, 5, 6, 7), (1, 8, 3, 9, 5, 10, 7, 11)])
    face = {1, 3, 5, 7}
    tris = []
    for a, b in m.source_cell_ranges:
        fs = []
        for t in m.tetrahedra[a:b]:
            for tri in __import__("itertools").combinations(t, 3):
                if set(tri) <= face:
                    fs.append(frozenset(tri))
        tris.append(set(fs))
    assert tris[0] == tris[1] and len(tris[0]) == 2


def test_p6_exhaustive_realizable_cube_face_diagonal_configurations():
    import itertools

    faces = ((0, 1, 3, 2), (4, 5, 7, 6), (0, 1, 5, 4), (2, 3, 7, 6), (0, 2, 6, 4), (1, 3, 7, 5))
    opposites = {frozenset(face): {face[i]: face[(i + 2) % 4] for i in range(4)} for face in faces}
    representatives = {}
    for ids in itertools.permutations(range(8)):
        pattern = []
        for face in faces:
            m = min(face, key=lambda v: ids[v])
            pattern.append(tuple(sorted((m, opposites[frozenset(face)][m]))))
        representatives.setdefault(tuple(pattern), ids)
    assert representatives
    coordinates = cube()
    for pattern, ids in representatives.items():
        result = tetrahedralize_cell(coordinates, vertex_ids=ids)
        assert len(result.tetrahedra) in (5, 6)
        counts = {}
        for tet in result.tetrahedra:
            for tri in itertools.combinations(tet, 3):
                key = frozenset(tri)
                counts[key] = counts.get(key, 0) + 1
        boundary = {tri for tri, count in counts.items() if count == 1}
        for face, expected in zip(faces, pattern, strict=True):
            global_face = {ids[v] for v in face}
            triangles = [tri for tri in boundary if tri <= global_face]
            assert len(triangles) == 2
            diagonal = tuple(sorted(set.intersection(*(set(t) for t in triangles))))
            assert diagonal == tuple(sorted(ids[v] for v in expected))


def test_p3_holes_disconnected_and_signed_winding_rules():
    outer = [(0, 0), (4, 0), (4, 4), (0, 4)]
    hole = [(1, 1), (1, 3), (3, 3), (3, 1)]
    r = polygonal_region_from_paths([outer, hole], fill_rule="nonzero")
    assert len(r.components) == 1 and len(r.components[0].hole_boundaries) == 1
    other = [(6, 0), (7, 0), (7, 1), (6, 1)]
    assert len(polygonal_region_from_paths([outer, other], fill_rule="positive").components) == 2
    assert (
        len(polygonal_region_from_paths([list(reversed(outer))], fill_rule="negative").components)
        == 1
    )


def test_p3_exact_algebraic_crossing():
    a = sp.sqrt(2)
    r = polygonal_region_from_paths([[(0, 0), (a, a), (0, a), (a, 0)]])
    vertices = {v for c in r.components for v in c.outer_boundary.vertices}
    assert (a / 2, a / 2) in vertices


def test_p4_convex_polytope_boundary_reconstruction():
    boundary = convert_region(Polytope(cube()), "boundary")
    assert isinstance(boundary, Polyhedron)
    assert len(boundary.components) == 1
    assert len(boundary.components[0].outer_shell.faces) == 6


def test_p2_general_zonotope_recognition_when_replay_certifies():
    from semialg.standard_regions import Zonotope

    generators = ((1, 0), (0, 1), (1, 1))
    points = []
    import itertools

    for bits in itertools.product((0, 1), repeat=3):
        points.append(tuple(sum(bits[k] * generators[k][j] for k in range(3)) for j in range(2)))
    result = canonicalize_polyhedron(Polytope(points))
    assert isinstance(result, Zonotope)


def test_p2_shell_nesting_reclassified_as_cavity():
    def shell(lo, hi):
        v = [
            (lo, lo, lo),
            (hi, lo, lo),
            (lo, hi, lo),
            (hi, hi, lo),
            (lo, lo, hi),
            (hi, lo, hi),
            (lo, hi, hi),
            (hi, hi, hi),
        ]
        f = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
        return PolyhedralShell(v, f)

    raw = Polyhedron((PolyhedralComponent(shell(0, 4)), PolyhedralComponent(shell(1, 2))))
    result = canonicalize_polyhedron(raw)
    assert isinstance(result, Polyhedron) and len(result.components) == 1
    assert len(result.components[0].cavity_shells) == 1
