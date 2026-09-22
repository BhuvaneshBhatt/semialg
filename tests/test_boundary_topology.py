import pytest

from semialg.boundary_topology import (
    PolygonalComponent,
    PolygonalSet,
    PolyhedralComponent,
    PolyhedralShell,
    Polyhedron,
    deduplicate_indexed_vertices,
    inner_polygons,
    inner_polyhedra,
    outer_polygons,
    outer_polyhedra,
    polygon_vertices,
    polyhedron_face_indices,
    polyhedron_vertices,
)
from semialg.standard_regions import Hexahedron, Polygon, Polytope

CUBE_VERTICES = (
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
)
CUBE_FACES = (
    (0, 3, 2, 1),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
    (4, 5, 6, 7),
)


def test_polygonal_set_orients_outer_and_hole_boundaries():
    component = PolygonalComponent(
        ((0, 0), (0, 4), (4, 4), (4, 0)),
        (((1, 1), (2, 1), (2, 2), (1, 2)),),
    )
    region = PolygonalSet((component, Polygon(((10, 0), (12, 0), (11, 1)))))

    assert len(outer_polygons(region)) == 2
    assert len(inner_polygons(region)) == 1
    assert component.outer_boundary.vertices[0] == (0, 0)
    assert component.hole_boundaries[0].vertices[0] == (1, 1)
    assert region.dimension() == 2
    assert region.ambient_dimension() == 2
    assert len(polygon_vertices(region)) == 11


def test_exact_vertex_deduplication_reindexes_cells_without_tolerance():
    result = deduplicate_indexed_vertices(
        ((0, 0, 0), (1, 0, 0), (0, 0, 0), (0, 1, 0)),
        ((0, 1, 3), (2, 3, 1)),
    )
    assert result.vertices == ((0, 0, 0), (1, 0, 0), (0, 1, 0))
    assert result.cells == ((0, 1, 2), (0, 2, 1))
    assert result.original_to_unique == (0, 1, 0, 2)


def test_polyhedral_shell_requires_closed_consistently_oriented_manifold():
    shell = PolyhedralShell(CUBE_VERTICES, CUBE_FACES)
    assert len(shell.edge_indices) == 12
    assert polyhedron_face_indices(shell) == (CUBE_FACES,)

    with pytest.raises(ValueError, match="closed 2-manifold"):
        PolyhedralShell(CUBE_VERTICES, CUBE_FACES[:-1])

    flipped = list(CUBE_FACES)
    flipped[1] = tuple(reversed(flipped[1]))
    with pytest.raises(ValueError, match="consistent orientation"):
        PolyhedralShell(CUBE_VERTICES, flipped)


def test_polyhedron_keeps_components_and_cavities_explicit():
    outer = PolyhedralShell(CUBE_VERTICES, CUBE_FACES)
    inner_vertices = tuple(tuple(value / 2 + 1 / 4 for value in vertex) for vertex in CUBE_VERTICES)
    # A cavity has the opposite boundary orientation from a filled outer shell.
    inner = PolyhedralShell(inner_vertices, tuple(tuple(reversed(face)) for face in CUBE_FACES))
    solid = Polyhedron((PolyhedralComponent(outer, (inner,)),))

    assert outer_polyhedra(solid)[0].signed_volume > 0
    assert inner_polyhedra(solid)[0].signed_volume < 0
    assert solid.dimension() == 3
    assert solid.ambient_dimension() == 3
    assert len(polyhedron_vertices(solid)) == 16


def test_hexahedron_is_a_canonical_polytope_type():
    hexahedron = Hexahedron(CUBE_VERTICES)
    assert isinstance(hexahedron, Polytope)
    assert len(hexahedron.vertices) == 8
    assert len(hexahedron.face_vertex_indices) == 6
    assert all(len(face) == 4 for face in hexahedron.face_vertex_indices)
