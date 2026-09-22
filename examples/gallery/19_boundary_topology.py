"""Canonical boundary topology: polygon holes and polyhedral cavities."""

from semialg import (
    PolygonalComponent,
    PolygonalSet,
    PolyhedralComponent,
    PolyhedralShell,
    Polyhedron,
    inner_polygons,
    inner_polyhedra,
    outer_polygons,
    outer_polyhedra,
)

polygonal_set = PolygonalSet(
    (
        PolygonalComponent(
            ((0, 0), (4, 0), (4, 4), (0, 4)),
            (((1, 1), (1, 2), (2, 2), (2, 1)),),
        ),
    )
)
print("outer polygon boundaries:", outer_polygons(polygonal_set))
print("polygon holes:", inner_polygons(polygonal_set))

cube_vertices = (
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
)
cube_faces = (
    (0, 3, 2, 1),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
    (4, 5, 6, 7),
)
outer_shell = PolyhedralShell(cube_vertices, cube_faces)
polyhedron = Polyhedron((PolyhedralComponent(outer_shell),))
print("outer polyhedral shells:", outer_polyhedra(polyhedron))
print("cavity shells:", inner_polyhedra(polyhedron))
