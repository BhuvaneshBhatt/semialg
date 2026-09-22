"""Canonical boundary-topology objects for polygonal and polyhedral sets.

This module separates boundary incidence from convexity and triangulation.  A
``Polygon`` remains one simple planar polygon.  ``PolygonalSet`` records one or
more connected filled components, each with an outer boundary and zero or more
hole boundaries.  ``Polyhedron`` analogously records closed oriented boundary
shells, including cavity shells, without requiring a tetrahedralization.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import sympy as sp

from ._standard_region_geometry import _PointData, _signed_polygon_orientation, _sympify_point
from ._standard_regions_base import Geometry
from ._standard_regions_polyhedral import Polygon
from ._zero_testing import certified_equal, certified_zero
from .exact_arithmetic import compare_exact_reals


def _canonical_cycle_start(indices: tuple[int, ...]) -> tuple[int, ...]:
    """Rotate a nonempty cycle to its lexicographically least rotation."""

    rotations = tuple(indices[offset:] + indices[:offset] for offset in range(len(indices)))
    return min(rotations)


def _normalize_polygon_orientation(polygon: Polygon, *, counterclockwise: bool) -> Polygon:
    orientation = _signed_polygon_orientation(polygon.vertices)
    expected = 1 if counterclockwise else -1
    vertices = polygon.vertices if orientation == expected else tuple(reversed(polygon.vertices))
    # Make equality independent of the arbitrary first vertex while retaining orientation.
    start = min(
        range(len(vertices)),
        key=lambda index: tuple(sp.default_sort_key(value) for value in vertices[index]),
    )
    vertices = vertices[start:] + vertices[:start]
    return Polygon(vertices, construction_conditions=polygon.construction_conditions)


@dataclass(frozen=True)
class PolygonalComponent:
    """One connected polygonal component with explicitly represented holes.

    The outer boundary is stored counterclockwise and hole boundaries clockwise.
    All boundaries are simple :class:`~semialg.standard_regions.Polygon` objects.
    Spatial containment/disjointness certification belongs to polygon
    canonicalization; this class records already-classified
    boundary topology without guessing from coordinates.
    """

    outer_boundary: Polygon
    hole_boundaries: tuple[Polygon, ...] = ()

    def __init__(
        self,
        outer_boundary: Polygon | Sequence[Sequence[object]],
        hole_boundaries: Sequence[Polygon | Sequence[Sequence[object]]] = (),
    ):
        outer = outer_boundary if isinstance(outer_boundary, Polygon) else Polygon(outer_boundary)
        holes = tuple(
            hole if isinstance(hole, Polygon) else Polygon(hole) for hole in hole_boundaries
        )
        object.__setattr__(
            self, "outer_boundary", _normalize_polygon_orientation(outer, counterclockwise=True)
        )
        object.__setattr__(
            self,
            "hole_boundaries",
            tuple(_normalize_polygon_orientation(hole, counterclockwise=False) for hole in holes),
        )


@dataclass(frozen=True)
class PolygonalSet(Geometry):
    """Canonical boundary topology for a planar polygonal set.

    A set may have multiple disconnected filled components; each component may
    contain any number of explicitly classified holes.
    """

    components: tuple[PolygonalComponent, ...]

    def __init__(
        self,
        components: Sequence[
            PolygonalComponent
            | Polygon
            | Sequence[Sequence[object]]
            | tuple[
                Polygon | Sequence[Sequence[object]], Sequence[Polygon | Sequence[Sequence[object]]]
            ]
        ],
    ):
        normalized: list[PolygonalComponent] = []
        for component in components:
            if isinstance(component, PolygonalComponent):
                normalized.append(component)
            elif isinstance(component, Polygon):
                normalized.append(PolygonalComponent(component))
            elif (
                isinstance(component, tuple)
                and len(component) == 2
                and isinstance(component[1], Sequence)
            ):
                normalized.append(PolygonalComponent(component[0], component[1]))
            else:
                normalized.append(PolygonalComponent(component))
        object.__setattr__(self, "components", tuple(normalized))

    def dimension(self) -> int:
        return 2 if self.components else -1

    def ambient_dimension(self) -> int:
        return 2


@dataclass(frozen=True)
class IndexedVertexData:
    """Result of exact vertex deduplication and cell reindexing."""

    vertices: tuple[_PointData, ...]
    cells: tuple[tuple[int, ...], ...]
    original_to_unique: tuple[int, ...]


def deduplicate_indexed_vertices(
    vertices: Sequence[Sequence[object]],
    cells: Sequence[Sequence[int]],
) -> IndexedVertexData:
    """Deduplicate exact vertices and remap zero-based cell indices.

    Equality is SymPy structural equality after sympification; no numerical
    tolerance is introduced.  The first occurrence of each coordinate is kept.
    """

    points = tuple(_sympify_point(vertex) for vertex in vertices)
    if points and len({len(point) for point in points}) != 1:
        raise ValueError("all vertices must have the same ambient dimension")

    unique_vertices: list[_PointData] = []
    original_to_unique: list[int] = []
    for point in points:
        index = None
        for candidate_index, candidate in enumerate(unique_vertices):
            decisions = tuple(
                certified_equal(left, right) for left, right in zip(point, candidate, strict=True)
            )
            if all(decision is True for decision in decisions):
                index = candidate_index
                break
        if index is None:
            index = len(unique_vertices)
            unique_vertices.append(point)
        original_to_unique.append(index)

    remapped_cells: list[tuple[int, ...]] = []
    for cell in cells:
        indices = tuple(int(index) for index in cell)
        if any(index < 0 or index >= len(points) for index in indices):
            raise IndexError("cell index is outside the vertex array")
        remapped_cells.append(tuple(original_to_unique[index] for index in indices))

    return IndexedVertexData(
        tuple(unique_vertices), tuple(remapped_cells), tuple(original_to_unique)
    )


def _normalize_face(face: Sequence[int], vertex_count: int) -> tuple[int, ...]:
    indices = tuple(int(index) for index in face)
    if len(indices) > 1 and indices[0] == indices[-1]:
        indices = indices[:-1]
    if len(indices) < 3:
        raise ValueError("a polyhedral face requires at least three vertices")
    if len(set(indices)) != len(indices):
        raise ValueError("a polyhedral face cannot repeat a vertex")
    if any(index < 0 or index >= vertex_count for index in indices):
        raise IndexError("face index is outside the vertex array")
    return indices


def _face_edges(face: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple((face[index], face[(index + 1) % len(face)]) for index in range(len(face)))


def _validate_closed_oriented_manifold(faces: tuple[tuple[int, ...], ...]) -> None:
    if not faces:
        raise ValueError("a polyhedral shell requires at least one face")

    undirected_counts: Counter[tuple[int, int]] = Counter()
    directed_counts: Counter[tuple[int, int]] = Counter()
    incident_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    for face_index, face in enumerate(faces):
        for start, end in _face_edges(face):
            edge = (min(start, end), max(start, end))
            undirected_counts[edge] += 1
            directed_counts[(start, end)] += 1
            incident_faces[edge].append(face_index)

    bad_edges = tuple(edge for edge, count in undirected_counts.items() if count != 2)
    if bad_edges:
        raise ValueError(
            "polyhedral shell must be a closed 2-manifold; every edge must have exactly two incident faces"
        )
    for start, end in undirected_counts:
        if directed_counts[(start, end)] != 1 or directed_counts[(end, start)] != 1:
            raise ValueError(
                "polyhedral shell faces must have consistent orientation across shared edges"
            )

    adjacency: dict[int, set[int]] = defaultdict(set)
    for adjacent in incident_faces.values():
        first, second = adjacent
        adjacency[first].add(second)
        adjacency[second].add(first)
    visited = {0}
    queue = deque([0])
    while queue:
        face_index = queue.popleft()
        for neighbor in adjacency[face_index]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    if len(visited) != len(faces):
        raise ValueError("a PolyhedralShell must be connected; use multiple shells for components")


def _newell_normal(vertices: tuple[_PointData, ...], face: tuple[int, ...]) -> sp.Matrix:
    normal = sp.zeros(3, 1)
    for index, start_index in enumerate(face):
        end_index = face[(index + 1) % len(face)]
        start = vertices[start_index]
        end = vertices[end_index]
        normal[0] += (start[1] - end[1]) * (start[2] + end[2])
        normal[1] += (start[2] - end[2]) * (start[0] + end[0])
        normal[2] += (start[0] - end[0]) * (start[1] + end[1])
    return normal


def _validate_planar_faces(
    vertices: tuple[_PointData, ...], faces: tuple[tuple[int, ...], ...]
) -> None:
    for face in faces:
        anchor = sp.Matrix(vertices[face[0]])
        directions = [sp.Matrix(vertices[index]) - anchor for index in face[1:]]
        if sp.Matrix.hstack(*directions).rank() != 2:
            raise ValueError("each polyhedral face must span an affine plane")
        normal = _newell_normal(vertices, face)
        if normal == sp.zeros(3, 1):
            raise ValueError("polyhedral face orientation is degenerate")
        for index in face[1:]:
            if certified_zero(normal.dot(sp.Matrix(vertices[index]) - anchor)) is not True:
                raise ValueError("polyhedral face vertices must be coplanar")


@dataclass(frozen=True)
class PolyhedralShell(Geometry):
    """One connected, closed, consistently oriented polygonal boundary shell."""

    vertices: tuple[_PointData, ...]
    faces: tuple[tuple[int, ...], ...]

    def __init__(
        self,
        vertices: Sequence[Sequence[object]],
        faces: Sequence[Sequence[int]],
        *,
        deduplicate_vertices: bool = True,
    ):
        raw_vertices = tuple(_sympify_point(vertex) for vertex in vertices)
        if not raw_vertices:
            raise ValueError("a polyhedral shell requires vertices")
        if any(len(vertex) != 3 for vertex in raw_vertices):
            raise ValueError("PolyhedralShell currently represents three-dimensional boundaries")

        raw_faces = tuple(tuple(int(index) for index in face) for face in faces)
        if deduplicate_vertices:
            deduplicated = deduplicate_indexed_vertices(raw_vertices, raw_faces)
            normalized_vertices = deduplicated.vertices
            candidate_faces = deduplicated.cells
        else:
            normalized_vertices = raw_vertices
            candidate_faces = raw_faces
        normalized_faces = tuple(
            _normalize_face(face, len(normalized_vertices)) for face in candidate_faces
        )
        if len(set(normalized_faces)) != len(normalized_faces):
            raise ValueError("a polyhedral shell cannot contain duplicate oriented faces")
        _validate_closed_oriented_manifold(normalized_faces)
        _validate_planar_faces(normalized_vertices, normalized_faces)
        object.__setattr__(self, "vertices", normalized_vertices)
        object.__setattr__(self, "faces", normalized_faces)

    def dimension(self) -> int:
        return 2

    def ambient_dimension(self) -> int:
        return 3

    @property
    def edge_indices(self) -> tuple[tuple[int, int], ...]:
        edges = {tuple(sorted(edge)) for face in self.faces for edge in _face_edges(face)}
        return tuple(sorted(edges))

    @property
    def signed_volume(self) -> sp.Expr:
        """Return the oriented volume enclosed by this closed shell."""

        six_volume = sp.Integer(0)
        for face in self.faces:
            anchor = sp.Matrix(self.vertices[face[0]])
            for offset in range(1, len(face) - 1):
                second = sp.Matrix(self.vertices[face[offset]])
                third = sp.Matrix(self.vertices[face[offset + 1]])
                six_volume += anchor.dot(second.cross(third))
        return sp.simplify(six_volume / 6)

    def with_reversed_orientation(self) -> PolyhedralShell:
        """Return the same shell with every face orientation reversed."""

        return PolyhedralShell(
            self.vertices,
            tuple(tuple(reversed(face)) for face in self.faces),
            deduplicate_vertices=False,
        )


@dataclass(frozen=True)
class PolyhedralComponent:
    """Boundary topology of one connected three-dimensional solid component."""

    outer_shell: PolyhedralShell
    cavity_shells: tuple[PolyhedralShell, ...] = ()

    def __init__(
        self,
        outer_shell: PolyhedralShell,
        cavity_shells: Sequence[PolyhedralShell] = (),
    ):
        if not isinstance(outer_shell, PolyhedralShell):
            raise TypeError("outer_shell must be a PolyhedralShell")
        cavities = tuple(cavity_shells)
        if any(not isinstance(shell, PolyhedralShell) for shell in cavities):
            raise TypeError("cavity_shells must contain PolyhedralShell objects")

        try:
            outer_sign = compare_exact_reals(outer_shell.signed_volume, sp.Integer(0))
        except (TypeError, ValueError, NotImplementedError) as exc:
            raise ValueError("outer-shell orientation could not be certified exactly") from exc
        if outer_sign == 0:
            raise ValueError("outer shell must enclose nonzero oriented volume")
        normalized_outer = (
            outer_shell if outer_sign > 0 else outer_shell.with_reversed_orientation()
        )

        normalized_cavities = []
        for shell in cavities:
            try:
                cavity_sign = compare_exact_reals(shell.signed_volume, sp.Integer(0))
            except (TypeError, ValueError, NotImplementedError) as exc:
                raise ValueError("cavity-shell orientation could not be certified exactly") from exc
            if cavity_sign == 0:
                raise ValueError("cavity shell must enclose nonzero oriented volume")
            normalized_cavities.append(
                shell if cavity_sign < 0 else shell.with_reversed_orientation()
            )

        object.__setattr__(self, "outer_shell", normalized_outer)
        object.__setattr__(self, "cavity_shells", tuple(normalized_cavities))


@dataclass(frozen=True)
class Polyhedron(Geometry):
    """Boundary representation of one or more three-dimensional solid components.

    Shell nesting is explicit: callers provide each component's outer shell and
    cavity shells. Nesting is established by canonicalization rather than inferred
    from raw shell geometry here.
    """

    components: tuple[PolyhedralComponent, ...]

    def __init__(self, components: Sequence[PolyhedralComponent | PolyhedralShell]):
        normalized = tuple(
            component
            if isinstance(component, PolyhedralComponent)
            else PolyhedralComponent(component)
            for component in components
        )
        object.__setattr__(self, "components", normalized)

    def dimension(self) -> int:
        return 3 if self.components else -1

    def ambient_dimension(self) -> int:
        return 3


def polygon_vertices(region: Polygon | PolygonalSet) -> tuple[_PointData, ...]:
    """Return unique polygon vertices in deterministic first-occurrence order."""

    if isinstance(region, Polygon):
        boundaries: Iterable[Polygon] = (region,)
    elif isinstance(region, PolygonalSet):
        boundaries = (
            boundary
            for component in region.components
            for boundary in (component.outer_boundary, *component.hole_boundaries)
        )
    else:
        raise TypeError("region must be a Polygon or PolygonalSet")
    return tuple(dict.fromkeys(vertex for boundary in boundaries for vertex in boundary.vertices))


def outer_polygons(region: Polygon | PolygonalSet) -> tuple[Polygon, ...]:
    """Return all outer polygon boundaries."""

    if isinstance(region, Polygon):
        return (_normalize_polygon_orientation(region, counterclockwise=True),)
    if isinstance(region, PolygonalSet):
        return tuple(component.outer_boundary for component in region.components)
    raise TypeError("region must be a Polygon or PolygonalSet")


def inner_polygons(region: Polygon | PolygonalSet) -> tuple[Polygon, ...]:
    """Return all polygon hole boundaries."""

    if isinstance(region, Polygon):
        return ()
    if isinstance(region, PolygonalSet):
        return tuple(hole for component in region.components for hole in component.hole_boundaries)
    raise TypeError("region must be a Polygon or PolygonalSet")


def polyhedron_vertices(region: PolyhedralShell | Polyhedron) -> tuple[_PointData, ...]:
    """Return unique boundary vertices in deterministic first-occurrence order."""

    if isinstance(region, PolyhedralShell):
        shells = (region,)
    elif isinstance(region, Polyhedron):
        shells = tuple(
            shell
            for component in region.components
            for shell in (component.outer_shell, *component.cavity_shells)
        )
    else:
        raise TypeError("region must be a PolyhedralShell or Polyhedron")
    return tuple(dict.fromkeys(vertex for shell in shells for vertex in shell.vertices))


def polyhedron_face_indices(
    region: PolyhedralShell | Polyhedron,
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Return face indices grouped by boundary shell.

    Each shell keeps its own local vertex indexing.  This avoids silently
    changing topology merely to manufacture one global coordinate array.
    """

    if isinstance(region, PolyhedralShell):
        return (region.faces,)
    if isinstance(region, Polyhedron):
        return tuple(
            shell.faces
            for component in region.components
            for shell in (component.outer_shell, *component.cavity_shells)
        )
    raise TypeError("region must be a PolyhedralShell or Polyhedron")


def outer_polyhedra(region: PolyhedralShell | Polyhedron) -> tuple[PolyhedralShell, ...]:
    """Return all outer boundary shells."""

    if isinstance(region, PolyhedralShell):
        return (region,)
    if isinstance(region, Polyhedron):
        return tuple(component.outer_shell for component in region.components)
    raise TypeError("region must be a PolyhedralShell or Polyhedron")


def inner_polyhedra(region: PolyhedralShell | Polyhedron) -> tuple[PolyhedralShell, ...]:
    """Return all cavity boundary shells."""

    if isinstance(region, PolyhedralShell):
        return ()
    if isinstance(region, Polyhedron):
        return tuple(shell for component in region.components for shell in component.cavity_shells)
    raise TypeError("region must be a PolyhedralShell or Polyhedron")


__all__ = [
    "IndexedVertexData",
    "PolygonalComponent",
    "PolygonalSet",
    "PolyhedralShell",
    "PolyhedralComponent",
    "Polyhedron",
    "deduplicate_indexed_vertices",
    "polygon_vertices",
    "outer_polygons",
    "inner_polygons",
    "polyhedron_vertices",
    "polyhedron_face_indices",
    "outer_polyhedra",
    "inner_polyhedra",
]
