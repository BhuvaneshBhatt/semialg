"""Exact polyhedral representations, face lattices, and incidence data."""

from __future__ import annotations

import itertools
from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ._zero_testing import certified_zero
from .exact_arithmetic import compare_exact_reals


def _point(data: Sequence[object]) -> tuple[sp.Expr, ...]:
    return tuple(sp.sympify(x) for x in data)


def _sign(expr: sp.Expr) -> int | None:
    value = sp.simplify(expr)
    if value == 0:
        return 0
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    try:
        return compare_exact_reals(value, sp.Integer(0))
    except (TypeError, ValueError, NotImplementedError):
        return None


def _affine_dimension(points: Sequence[Sequence[object]]) -> int:
    if not points:
        return -1
    if len(points) == 1:
        return 0
    anchor = sp.Matrix(points[0])
    return sp.Matrix.hstack(*(sp.Matrix(p) - anchor for p in points[1:])).rank()


@dataclass(frozen=True)
class PolytopeFacet:
    """One supporting facet of a full-dimensional convex polytope."""

    vertex_indices: tuple[int, ...]
    normal: tuple[sp.Expr, ...]
    offset: sp.Expr

    def contains_vertex(self, index: int) -> bool:
        return index in self.vertex_indices


@dataclass(frozen=True)
class PolytopeFace:
    """One nonempty face in a polytope face lattice."""

    dimension: int
    vertex_indices: tuple[int, ...]
    facet_indices: tuple[int, ...]

    def contains_vertex(self, index: int) -> bool:
        return index in self.vertex_indices


@dataclass(frozen=True)
class PolytopeIncidence:
    """Exact vertex-facet incidence relation for a convex polytope."""

    facets: tuple[PolytopeFacet, ...]
    vertex_facets: tuple[tuple[int, ...], ...]

    @property
    def facet_count(self) -> int:
        return len(self.facets)


@dataclass(frozen=True)
class PolytopeAdjacency:
    """Exact vertex and facet adjacency relations."""

    vertex_neighbors: tuple[tuple[int, ...], ...]
    facet_neighbors: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class PolytopeFaceLattice:
    """Finite nonempty face lattice of a full-dimensional convex polytope."""

    dimension: int
    faces: tuple[PolytopeFace, ...]

    def faces_of_dimension(self, dimension: int) -> tuple[PolytopeFace, ...]:
        return tuple(face for face in self.faces if face.dimension == dimension)

    @property
    def f_vector(self) -> tuple[int, ...]:
        """Return conventional face counts ``(f0, ..., f_(d-1))``."""

        return tuple(len(self.faces_of_dimension(k)) for k in range(self.dimension))

    @property
    def extended_f_vector(self) -> tuple[int, ...]:
        """Return face counts including the empty and whole-polytope faces."""

        return (1, *self.f_vector, 1)

    @property
    def euler_characteristic(self) -> int:
        """Return the Euler characteristic of the closed polytope."""

        counts = (*self.f_vector, 1)
        return sum((-1) ** k * count for k, count in enumerate(counts))

    @property
    def boundary_euler_characteristic(self) -> int:
        """Return the Euler characteristic of the polytope boundary."""

        return sum((-1) ** k * count for k, count in enumerate(self.f_vector))


@dataclass(frozen=True)
class HConstraintRedundancyCertificate:
    """Replayable redundancy decision for one inequality of an H-representation."""

    constraint_index: int
    redundant: bool
    remaining_indices: tuple[int, ...]


@dataclass(frozen=True)
class HRepresentation:
    """Closed halfspace representation ``A x <= b`` of a polyhedron."""

    matrix: sp.ImmutableMatrix
    offsets: sp.ImmutableMatrix

    def __init__(self, matrix: Sequence[Sequence[object]], offsets: Sequence[object]):
        A = sp.ImmutableMatrix([[sp.sympify(x) for x in row] for row in matrix])
        b = sp.ImmutableMatrix([sp.sympify(x) for x in offsets])
        if A.rows == 0 or A.cols == 0:
            raise ValueError("an H-representation requires at least one nontrivial halfspace")
        if b.cols != 1 or b.rows != A.rows:
            raise ValueError("halfspace offsets must match the matrix row count")
        for row in range(A.rows):
            zero_entries = tuple(certified_zero(A[row, col]) for col in range(A.cols))
            if all(value is True for value in zero_entries):
                raise ValueError("halfspace normals must be nonzero")
            if not any(value is False for value in zero_entries):
                raise ValueError("could not certify that a halfspace normal is nonzero")
        object.__setattr__(self, "matrix", A)
        object.__setattr__(self, "offsets", b)

    @property
    def ambient_dimension(self) -> int:
        return self.matrix.cols

    def as_formula(self, variables: Sequence[sp.Symbol]) -> sp.Expr:
        variables = tuple(variables)
        if len(variables) != self.ambient_dimension:
            raise ValueError("variable count must match the H-representation dimension")
        x = sp.Matrix(variables)
        return sp.And(
            *(
                sp.Le((self.matrix.row(i) * x)[0], self.offsets[i, 0])
                for i in range(self.matrix.rows)
            )
        )

    def vertices(self) -> tuple[tuple[sp.Expr, ...], ...]:
        """Enumerate exact vertices from active sets of ambient dimension."""

        n = self.ambient_dimension
        vertices: list[tuple[sp.Expr, ...]] = []
        for active in itertools.combinations(range(self.matrix.rows), n):
            A = self.matrix[list(active), :]
            if certified_zero(A.det()) is True:
                continue
            b = self.offsets[list(active), :]
            sol = A.inv() * b
            point = tuple(sp.simplify(sol[i, 0]) for i in range(n))
            valid = True
            for row in range(self.matrix.rows):
                delta = sp.simplify((self.matrix.row(row) * sol)[0] - self.offsets[row, 0])
                sign = _sign(delta)
                if sign is None:
                    raise NotImplementedError(
                        "could not certify an H-representation vertex inequality"
                    )
                if sign > 0:
                    valid = False
                    break
            if valid and point not in vertices:
                vertices.append(point)
        return tuple(vertices)

    def is_bounded(self) -> bool:
        """Return whether the represented polyhedron is bounded."""

        from .reasoning_regions import region_bounded

        vars_ = sp.symbols(f"x0:{self.ambient_dimension}", real=True)
        return region_bounded(self.as_formula(vars_), vars_)

    def redundancy_certificates(self) -> tuple[HConstraintRedundancyCertificate, ...]:
        """Return exact replayable redundancy decisions for all inequalities."""

        from .decision import is_satisfiable

        variables = sp.symbols(f"x0:{self.ambient_dimension}", real=True)
        x = sp.Matrix(variables)
        certificates = []
        for index in range(self.matrix.rows):
            remaining = tuple(i for i in range(self.matrix.rows) if i != index)
            constraints = [
                sp.Le((self.matrix.row(i) * x)[0], self.offsets[i, 0]) for i in remaining
            ]
            violation = sp.Gt((self.matrix.row(index) * x)[0], self.offsets[index, 0])
            satisfiable = is_satisfiable(sp.And(*constraints, violation), variables)
            certificates.append(
                HConstraintRedundancyCertificate(index, not bool(satisfiable), remaining)
            )
        return tuple(certificates)

    def irredundant(self) -> HRepresentation:
        """Return an equivalent H-representation with redundant rows removed."""

        current = self
        index = current.matrix.rows - 1
        while index >= 0:
            certificate = current.redundancy_certificates()[index]
            if certificate.redundant and current.matrix.rows > 1:
                keep = [i for i in range(current.matrix.rows) if i != index]
                current = HRepresentation(
                    [[current.matrix[i, j] for j in range(current.matrix.cols)] for i in keep],
                    [current.offsets[i, 0] for i in keep],
                )
            index -= 1
        return current

    def to_polytope(self):
        """Convert a nonempty bounded full-dimensional H-polyhedron to a polytope."""

        if not self.is_bounded():
            raise ValueError("halfspaces do not define a bounded polytope")
        vertices = self.vertices()
        if not vertices:
            raise ValueError("halfspaces do not define a nonempty full-dimensional polytope")
        from .standard_regions import Polytope

        result = Polytope(vertices)
        if result.dimension() != self.ambient_dimension:
            raise ValueError("halfspaces do not define a full-dimensional polytope")
        object.__setattr__(result, "_h_representation", self.irredundant())
        return result


def verify_h_redundancy_certificate(
    representation: HRepresentation, certificate: HConstraintRedundancyCertificate
) -> bool:
    """Replay one H-constraint redundancy certificate against its representation."""

    index = certificate.constraint_index
    if index < 0 or index >= representation.matrix.rows:
        return False
    expected = tuple(i for i in range(representation.matrix.rows) if i != index)
    if certificate.remaining_indices != expected:
        return False
    return representation.redundancy_certificates()[index] == certificate


def h_representation_from_vertices(
    vertices: Sequence[Sequence[object]],
) -> HRepresentation:
    """Compute the exact irredundant supporting halfspaces of a full-dimensional hull."""

    verts = tuple(_point(v) for v in vertices)
    if not verts:
        raise ValueError("a polytope requires vertices")
    n = len(verts[0])
    if any(len(v) != n for v in verts):
        raise ValueError("polytope vertices must have one ambient dimension")
    anchor = sp.Matrix(verts[0])
    if sp.Matrix.hstack(*(sp.Matrix(v) - anchor for v in verts[1:])).rank() != n:
        raise ValueError("H-representation currently requires a full-dimensional polytope")

    facets: dict[frozenset[int], tuple[tuple[sp.Expr, ...], sp.Expr]] = {}
    for combo in itertools.combinations(range(len(verts)), n):
        base = sp.Matrix(verts[combo[0]])
        differences = sp.Matrix.vstack(*(sp.Matrix(verts[i]).T - base.T for i in combo[1:]))
        null = differences.nullspace()
        if len(null) != 1:
            continue
        normal = null[0]
        offset = sp.simplify((normal.T * base)[0])
        deltas = [sp.simplify((normal.T * sp.Matrix(v))[0] - offset) for v in verts]
        signs = [_sign(delta) for delta in deltas]
        if any(sign is None for sign in signs):
            raise NotImplementedError("could not certify a supporting hyperplane")
        nonzero = {sign for sign in signs if sign != 0}
        if len(nonzero) > 1 or not nonzero:
            continue
        if 1 in nonzero:
            normal = -normal
            offset = -offset
            deltas = [-d for d in deltas]
        incident = frozenset(i for i, d in enumerate(deltas) if certified_zero(d) is True)
        incident_points = [sp.Matrix(verts[i]) for i in incident]
        facet_anchor = incident_points[0]
        rank = sp.Matrix.hstack(*(p - facet_anchor for p in incident_points[1:])).rank()
        if rank != n - 1:
            continue
        facets[incident] = (tuple(sp.simplify(x) for x in normal), sp.simplify(offset))

    if not facets:
        raise ValueError("could not derive polytope facets")
    ordered = sorted(facets.items(), key=lambda item: tuple(sorted(item[0])))
    return HRepresentation(
        [normal for _, (normal, _) in ordered],
        [offset for _, (_, offset) in ordered],
    )


def polytope_incidence(vertices: Sequence[Sequence[object]]) -> PolytopeIncidence:
    """Return exact facet and vertex-facet incidence for a full-dimensional polytope."""

    verts = tuple(_point(v) for v in vertices)
    hrep = h_representation_from_vertices(verts)
    facets: list[PolytopeFacet] = []
    for row in range(hrep.matrix.rows):
        normal = tuple(hrep.matrix[row, col] for col in range(hrep.matrix.cols))
        offset = hrep.offsets[row, 0]
        incident = tuple(
            i
            for i, vertex in enumerate(verts)
            if certified_zero(sum(a * x for a, x in zip(normal, vertex, strict=True)) - offset)
            is True
        )
        facets.append(PolytopeFacet(incident, normal, offset))
    vertex_facets = tuple(
        tuple(j for j, facet in enumerate(facets) if i in facet.vertex_indices)
        for i in range(len(verts))
    )
    return PolytopeIncidence(tuple(facets), vertex_facets)


def polytope_face_lattice(vertices: Sequence[Sequence[object]]) -> PolytopeFaceLattice:
    """Construct the complete nonempty face lattice from exact facet incidence."""

    verts = tuple(_point(v) for v in vertices)
    incidence = polytope_incidence(verts)
    all_vertices = frozenset(range(len(verts)))
    face_sets: set[frozenset[int]] = {all_vertices}
    face_sets.update(frozenset(f.vertex_indices) for f in incidence.facets)
    changed = True
    while changed:
        changed = False
        current = tuple(face_sets)
        for left, right in itertools.combinations(current, 2):
            intersection = left & right
            if intersection and intersection not in face_sets:
                face_sets.add(intersection)
                changed = True

    faces = []
    facet_sets = tuple(frozenset(f.vertex_indices) for f in incidence.facets)
    for vertex_set in face_sets:
        points = tuple(verts[i] for i in sorted(vertex_set))
        dimension = _affine_dimension(points)
        containing = tuple(i for i, facet in enumerate(facet_sets) if vertex_set <= facet)
        faces.append(PolytopeFace(dimension, tuple(sorted(vertex_set)), containing))
    faces.sort(key=lambda face: (face.dimension, face.vertex_indices))
    dimension = _affine_dimension(verts)
    return PolytopeFaceLattice(dimension, tuple(faces))


def polytope_adjacency(vertices: Sequence[Sequence[object]]) -> PolytopeAdjacency:
    """Return exact vertex-edge and facet-ridge adjacency."""

    lattice = polytope_face_lattice(vertices)
    edges = lattice.faces_of_dimension(1)
    vertex_neighbors = [set() for _ in vertices]
    for edge in edges:
        if len(edge.vertex_indices) < 2:
            continue
        for i, j in itertools.combinations(edge.vertex_indices, 2):
            vertex_neighbors[i].add(j)
            vertex_neighbors[j].add(i)

    facet_count = len(lattice.faces_of_dimension(lattice.dimension - 1))
    facet_neighbors = [set() for _ in range(facet_count)]
    for ridge in lattice.faces_of_dimension(lattice.dimension - 2):
        for i, j in itertools.combinations(ridge.facet_indices, 2):
            facet_neighbors[i].add(j)
            facet_neighbors[j].add(i)
    return PolytopeAdjacency(
        tuple(tuple(sorted(x)) for x in vertex_neighbors),
        tuple(tuple(sorted(x)) for x in facet_neighbors),
    )


def polytopes_combinatorially_equivalent(left, right) -> bool:
    """Return whether two full-dimensional polytopes have isomorphic face lattices."""

    if left.dimension() != right.dimension():
        return False
    if len(left.vertices) != len(right.vertices):
        return False
    left_lattice = polytope_face_lattice(left.vertices)
    right_lattice = polytope_face_lattice(right.vertices)
    if left_lattice.f_vector != right_lattice.f_vector:
        return False

    left_inc = polytope_incidence(left.vertices)
    right_inc = polytope_incidence(right.vertices)
    if len(left_inc.facets) != len(right_inc.facets):
        return False
    target_facets = {frozenset(f.vertex_indices) for f in right_inc.facets}
    left_adj = polytope_adjacency(left.vertices).vertex_neighbors
    right_adj = polytope_adjacency(right.vertices).vertex_neighbors

    left_sig = [
        (len(left_adj[i]), tuple(sorted(len(left_inc.facets[j].vertex_indices) for j in fs)))
        for i, fs in enumerate(left_inc.vertex_facets)
    ]
    right_sig = [
        (len(right_adj[i]), tuple(sorted(len(right_inc.facets[j].vertex_indices) for j in fs)))
        for i, fs in enumerate(right_inc.vertex_facets)
    ]
    candidates = [
        tuple(j for j, sig in enumerate(right_sig) if sig == left_sig[i])
        for i in range(len(left.vertices))
    ]
    if any(not group for group in candidates):
        return False

    order = sorted(range(len(left.vertices)), key=lambda i: len(candidates[i]))
    mapping: dict[int, int] = {}
    used: set[int] = set()

    def search(position: int) -> bool:
        if position == len(order):
            mapped = {
                frozenset(mapping[i] for i in facet.vertex_indices) for facet in left_inc.facets
            }
            return mapped == target_facets
        source = order[position]
        for target in candidates[source]:
            if target in used:
                continue
            if any(
                ((other in left_adj[source]) != (mapping[other] in right_adj[target]))
                for other in mapping
            ):
                continue
            mapping[source] = target
            used.add(target)
            if search(position + 1):
                return True
            used.remove(target)
            del mapping[source]
        return False

    return search(0)


__all__ = [
    "HConstraintRedundancyCertificate",
    "HRepresentation",
    "PolytopeAdjacency",
    "PolytopeFace",
    "PolytopeFaceLattice",
    "PolytopeFacet",
    "PolytopeIncidence",
    "h_representation_from_vertices",
    "polytope_adjacency",
    "polytope_face_lattice",
    "polytope_incidence",
    "polytopes_combinatorially_equivalent",
    "verify_h_redundancy_certificate",
]
