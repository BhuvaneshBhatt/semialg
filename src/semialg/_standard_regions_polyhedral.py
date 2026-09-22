# ruff: noqa: F403, F405
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import sympy as sp

from ._standard_region_geometry import (
    _affine_dimension,
    _independent_vectors,
    _PointData,
    _provably_nonpositive,
    _provably_zero,
    _signed_polygon_orientation,
    _sympify_point,
    _triangulate_polygon,
    _validate_interval,
    _validate_same_dimension,
    _validate_simple_polygon,
    _validate_triangle_angle,
    _validate_triangle_sides,
    _vector_rank,
)
from ._standard_regions_base import *
from ._standard_regions_base import _validate_vector


@dataclass(frozen=True)
class Interval(StandardRegion):
    lower: sp.Expr
    upper: sp.Expr
    lower_closed: bool = True
    upper_closed: bool = True

    def __init__(
        self, lower: object, upper: object, *, lower_closed: bool = True, upper_closed: bool = True
    ):
        lower_expr = sp.sympify(lower)
        upper_expr = sp.sympify(upper)
        _validate_interval(lower_expr, upper_expr, label="interval bounds")
        object.__setattr__(self, "lower", lower_expr)
        object.__setattr__(self, "upper", upper_expr)
        object.__setattr__(self, "lower_closed", lower_closed)
        object.__setattr__(self, "upper_closed", upper_closed)

    def dimension(self) -> int:
        if not _provably_zero(self.upper - self.lower):
            return 1
        return 0 if self.lower_closed and self.upper_closed else -1

    def ambient_dimension(self) -> int:
        return 1


@dataclass(frozen=True)
class Box(StandardRegion):
    bounds: tuple[tuple[sp.Expr, sp.Expr], ...]

    def __init__(self, bounds: Sequence[tuple[object, object]]):
        normalized = tuple((sp.sympify(a), sp.sympify(b)) for a, b in bounds)
        for lower, upper in normalized:
            _validate_interval(lower, upper, label="box bounds")
        object.__setattr__(self, "bounds", normalized)

    def dimension(self) -> int:
        return sum(not _provably_zero(upper - lower) for lower, upper in self.bounds)

    def ambient_dimension(self) -> int:
        return len(self.bounds)


@dataclass(frozen=True)
class _SimplexBase(StandardRegion):
    vertices: tuple[_PointData, ...]

    def __init__(self, vertices: Sequence[Sequence[object]]):
        verts = tuple(_sympify_point(v) for v in vertices)
        if not verts:
            raise ValueError("a simplex requires at least one vertex")
        _validate_same_dimension(verts, label="simplex vertices")
        object.__setattr__(self, "vertices", verts)

    def dimension(self) -> int:
        return _affine_dimension(self.vertices)

    def ambient_dimension(self) -> int:
        return len(self.vertices[0]) if self.vertices else 0


@dataclass(frozen=True)
class _PolygonBase(StandardRegion):
    vertices: tuple[_PointData, ...]

    def __init__(self, vertices: Sequence[Sequence[object]]):
        verts = tuple(_sympify_point(v) for v in vertices)
        if len(verts) > 1 and verts[0] == verts[-1]:
            verts = verts[:-1]
        if len(verts) < 3:
            raise ValueError("a polygon requires at least three vertices")
        if len({len(v) for v in verts}) != 1 or len(verts[0]) != 2:
            raise ValueError("_PolygonBase represents 2D polygons")
        _validate_simple_polygon(verts)
        _signed_polygon_orientation(verts)
        object.__setattr__(self, "vertices", verts)

    def dimension(self) -> int:
        return 2

    def ambient_dimension(self) -> int:
        return 2

    def triangulation(self) -> tuple[Simplex, ...]:
        return tuple(Simplex(triangle) for triangle in _triangulate_polygon(self.vertices))


@dataclass(frozen=True)
class Polytope(StandardRegion):
    """Convex polytope represented canonically by its extreme vertices."""

    vertices: tuple[_PointData, ...]
    _h_representation: object | None = field(default=None, compare=False, repr=False)

    def __init__(self, vertices: Sequence[Sequence[object]]):
        verts = tuple(_sympify_point(v) for v in vertices)
        if not verts:
            raise ValueError("a polytope requires at least one vertex")
        _validate_same_dimension(verts, label="polytope vertices")
        # Preserve first occurrence order while removing exact duplicate vertices.
        unique = tuple(dict.fromkeys(verts))
        object.__setattr__(self, "vertices", unique)
        object.__setattr__(self, "_h_representation", None)

    def dimension(self) -> int:
        return _affine_dimension(self.vertices)

    def ambient_dimension(self) -> int:
        return len(self.vertices[0])

    def affine_hull(self):
        anchor = self.vertices[0]
        directions = tuple(
            tuple(sp.Matrix(vertex) - sp.Matrix(anchor)) for vertex in self.vertices[1:]
        )
        return AffineSpace(anchor, directions)

    @classmethod
    def from_halfspaces(
        cls, matrix: Sequence[Sequence[object]], offsets: Sequence[object]
    ) -> Polytope:
        """Construct a bounded full-dimensional polytope from ``A x <= b``."""

        from .polyhedral import HRepresentation

        return HRepresentation(matrix, offsets).to_polytope()

    def h_representation(self):
        """Return the exact irredundant supporting-halfspace representation."""

        from .polyhedral import h_representation_from_vertices

        if self._h_representation is not None:
            return self._h_representation
        return h_representation_from_vertices(self.vertices)

    def incidence(self):
        """Return exact vertex-facet incidence data."""

        from .polyhedral import polytope_incidence

        return polytope_incidence(self.vertices)

    def facets(self):
        """Return the exact facets of this full-dimensional polytope."""

        return self.incidence().facets

    def face_lattice(self):
        """Return the complete exact nonempty face lattice."""

        from .polyhedral import polytope_face_lattice

        return polytope_face_lattice(self.vertices)

    def faces(self, dimension: int | None = None):
        """Return all nonempty faces, optionally restricted by dimension."""

        lattice = self.face_lattice()
        return lattice.faces if dimension is None else lattice.faces_of_dimension(dimension)

    def edges(self):
        """Return the one-dimensional faces of the polytope."""

        return self.faces(1)

    def ridges(self):
        """Return codimension-two faces of the polytope."""

        return self.faces(self.dimension() - 2)

    def adjacency(self):
        """Return exact vertex and facet adjacency relations."""

        from .polyhedral import polytope_adjacency

        return polytope_adjacency(self.vertices)

    def f_vector(self) -> tuple[int, ...]:
        """Return conventional face counts ``(f0, ..., f_(d-1))``."""

        return self.face_lattice().f_vector

    def euler_characteristic(self) -> int:
        """Return the Euler characteristic of this closed convex polytope."""

        return self.face_lattice().euler_characteristic

    def is_combinatorially_equivalent(self, other) -> bool:
        """Return whether ``other`` has an isomorphic face lattice."""

        if not isinstance(other, Polytope):
            return False
        from .polyhedral import polytopes_combinatorially_equivalent

        return polytopes_combinatorially_equivalent(self, other)

    def validate_topology(
        self,
        *,
        vertex_count: int | None = None,
        facet_count: int | None = None,
        facet_sizes: Sequence[int] | None = None,
    ) -> bool:
        """Validate requested exact hull topology invariants."""

        data = self.incidence()
        extreme = {i for facet in data.facets for i in facet.vertex_indices}
        if len(extreme) != len(self.vertices):
            return False
        if vertex_count is not None and len(self.vertices) != vertex_count:
            return False
        if facet_count is not None and data.facet_count != facet_count:
            return False
        if facet_sizes is not None and sorted(len(f.vertex_indices) for f in data.facets) != sorted(
            facet_sizes
        ):
            return False
        return True

    @classmethod
    def from_region(cls, region: StandardRegion) -> Polytope:
        """Convert a vertex-defined standard region to a canonical polytope."""

        if isinstance(region, (_SimplexBase, _PolygonBase)):
            return cls(region.vertices)
        if isinstance(region, _ParallelepipedBase):
            vertices = []
            for mask in __import__("itertools").product((0, 1), repeat=len(region.vectors)):
                point = tuple(
                    region.origin[j]
                    + sum(bit * vec[j] for bit, vec in zip(mask, region.vectors, strict=True))
                    for j in range(region.ambient_dimension())
                )
                vertices.append(point)
            return cls(vertices)
        raise TypeError(f"cannot convert {type(region).__name__} to Polytope")


@dataclass(frozen=True)
class Simplex(_SimplexBase):
    """Canonical simplex defined by affinely independent vertices."""

    construction_conditions: tuple[sp.Expr, ...] = field(default=(), compare=False)

    def __init__(
        self,
        vertices: Sequence[Sequence[object]] | _SimplexBase,
        *,
        construction_conditions: Sequence[object] = (),
    ):
        data = vertices.vertices if isinstance(vertices, _SimplexBase) else vertices
        super().__init__(data)
        conditions = tuple(
            condition
            for condition in (sp.simplify(sp.sympify(item)) for item in construction_conditions)
            if condition is not sp.true
        )
        object.__setattr__(self, "construction_conditions", conditions)

    def is_nondegenerate(self) -> bool:
        """Return whether the vertices are affinely independent."""

        return len(self.vertices) == self.dimension() + 1

    def affine_hull(self):
        anchor = self.vertices[0]
        directions = tuple(
            tuple(sp.Matrix(vertex) - sp.Matrix(anchor)) for vertex in self.vertices[1:]
        )
        return AffineSpace(anchor, directions)

    @classmethod
    def from_region(cls, region: _SimplexBase) -> Simplex:
        return cls(region)


class Triangle:
    """Convenience constructor namespace for canonical two-dimensional simplices.

    ``Triangle(...)`` and the named constructors return :class:`Simplex`
    instances; ``Triangle`` is not a separate runtime geometry type.
    """

    def __new__(cls, vertices: Sequence[Sequence[object]]) -> Simplex:
        return Simplex(vertices)

    @staticmethod
    def from_sides(a: object, b: object, c: object) -> Simplex:
        """Construct a triangle with side lengths ``a``, ``b``, and ``c``.

        The returned placement uses the side of length ``c`` on the positive
        x-axis, with lengths ``a`` and ``b`` opposite the first and second
        vertices respectively. Provably invalid symbolic data are rejected.
        """

        a, b, c = map(sp.sympify, (a, b, c))
        _validate_triangle_sides(a, b, c)
        x = sp.simplify((b**2 + c**2 - a**2) / (2 * c))
        y2 = sp.factor(b**2 - x**2)
        y = sp.sqrt(y2)
        return Simplex(
            ((0, 0), (c, 0), (x, y)),
            construction_conditions=(
                a > 0,
                b > 0,
                c > 0,
                a + b > c,
                a + c > b,
                b + c > a,
            ),
        )

    @staticmethod
    def from_sas(a: object, angle: object, b: object) -> Simplex:
        """Construct a triangle from two sides and their included angle."""

        a, angle, b = map(sp.sympify, (a, angle, b))
        if _provably_nonpositive(a) or _provably_nonpositive(b):
            raise ValueError("triangle side lengths must be positive")
        _validate_triangle_angle(angle)
        return Simplex(
            ((0, 0), (a, 0), (b * sp.cos(angle), b * sp.sin(angle))),
            construction_conditions=(a > 0, b > 0, angle > 0, angle < sp.pi),
        )

    @staticmethod
    def from_asa(angle_a: object, side: object, angle_b: object) -> Simplex:
        """Construct a triangle from two angles and their included side."""

        angle_a, side, angle_b = map(sp.sympify, (angle_a, side, angle_b))
        if _provably_nonpositive(side):
            raise ValueError("triangle side lengths must be positive")
        _validate_triangle_angle(angle_a)
        _validate_triangle_angle(angle_b)
        angle_c = sp.simplify(sp.pi - angle_a - angle_b)
        _validate_triangle_angle(angle_c)
        # Put the included side AB on the x-axis. The two rays meet at C.
        sin_c = sp.sin(angle_c)
        ac = sp.simplify(side * sp.sin(angle_b) / sin_c)
        return Simplex(
            ((0, 0), (side, 0), (ac * sp.cos(angle_a), ac * sp.sin(angle_a))),
            construction_conditions=(
                side > 0,
                angle_a > 0,
                angle_b > 0,
                angle_c > 0,
            ),
        )

    @staticmethod
    def from_aas(angle_a: object, angle_b: object, opposite_a: object) -> Simplex:
        """Construct a triangle from two angles and the side opposite the first."""

        angle_a, angle_b, opposite_a = map(sp.sympify, (angle_a, angle_b, opposite_a))
        if _provably_nonpositive(opposite_a):
            raise ValueError("triangle side lengths must be positive")
        _validate_triangle_angle(angle_a)
        _validate_triangle_angle(angle_b)
        angle_c = sp.simplify(sp.pi - angle_a - angle_b)
        _validate_triangle_angle(angle_c)
        scale = sp.simplify(opposite_a / sp.sin(angle_a))
        side_c = sp.simplify(scale * sp.sin(angle_c))
        side_b = sp.simplify(scale * sp.sin(angle_b))
        vertex = (side_b * sp.cos(angle_a), side_b * sp.sin(angle_a))
        return Simplex(
            ((0, 0), (side_c, 0), vertex),
            construction_conditions=(
                opposite_a > 0,
                angle_a > 0,
                angle_b > 0,
                angle_c > 0,
            ),
        )


def RegularPolygon(
    sides: int,
    center: Sequence[object] = (0, 0),
    radius: object = 1,
    *,
    rotation: object = 0,
) -> Polygon:
    """Return a canonical regular polygon as a :class:`Polygon`.

    Vertices lie on the circumcircle with the first vertex at ``rotation``
    radians from the positive x-axis.
    """

    if isinstance(sides, bool) or not isinstance(sides, int):
        raise TypeError("regular polygon side count must be an integer")
    if sides < 3:
        raise ValueError("a regular polygon requires at least three sides")
    center_pt = _sympify_point(center)
    if len(center_pt) != 2:
        raise ValueError("RegularPolygon requires a two-dimensional center")
    r = sp.sympify(radius)
    if _provably_nonpositive(r):
        raise ValueError("regular polygon radius must be positive")
    theta0 = sp.sympify(rotation)
    vertices = tuple(
        (
            sp.simplify(center_pt[0] + r * sp.cos(theta0 + 2 * sp.pi * k / sides)),
            sp.simplify(center_pt[1] + r * sp.sin(theta0 + 2 * sp.pi * k / sides)),
        )
        for k in range(sides)
    )
    return Polygon(vertices, construction_conditions=(r > 0,))


@dataclass(frozen=True)
class Polygon(_PolygonBase):
    """Canonical simple polygon in two-dimensional affine coordinates."""

    construction_conditions: tuple[sp.Expr, ...] = field(default=(), compare=False)

    def __init__(
        self,
        vertices: Sequence[Sequence[object]] | _PolygonBase,
        *,
        construction_conditions: Sequence[object] = (),
    ):
        data = vertices.vertices if isinstance(vertices, _PolygonBase) else vertices
        super().__init__(data)
        conditions = tuple(
            condition
            for condition in (sp.simplify(sp.sympify(item)) for item in construction_conditions)
            if condition is not sp.true
        )
        object.__setattr__(self, "construction_conditions", conditions)

    def affine_hull(self):
        return AffineSpace(self.vertices[0], ((1, 0), (0, 1)))

    @classmethod
    def from_region(cls, region: _PolygonBase) -> Polygon:
        return cls(region)


@dataclass(frozen=True)
class PolyhedralCone(StandardRegion):
    """Affine polyhedral cone with lineality directions and nonnegative generating rays."""

    point: _PointData
    directions: tuple[_PointData, ...]
    rays: tuple[_PointData, ...]

    def __init__(
        self,
        point: Sequence[object],
        directions: Sequence[Sequence[object]] = (),
        rays: Sequence[Sequence[object]] = (),
    ):
        anchor = _sympify_point(point)
        dirs = tuple(_sympify_point(v) for v in directions)
        ray_data = tuple(_validate_vector(v, label="polyhedral-cone ray") for v in rays)
        if not anchor:
            raise ValueError("a polyhedral cone requires a nonempty point")
        _validate_same_dimension((anchor, *dirs, *ray_data), label="polyhedral-cone data")
        basis = _independent_vectors(
            tuple(v for v in dirs if any(not _provably_zero(x) for x in v))
        )
        object.__setattr__(self, "point", anchor)
        object.__setattr__(self, "directions", basis)
        object.__setattr__(self, "rays", ray_data)

    def dimension(self) -> int:
        return _vector_rank((*self.directions, *self.rays))

    def ambient_dimension(self) -> int:
        return len(self.point)

    def affine_hull(self):
        return AffineSpace(self.point, (*self.directions, *self.rays))


@dataclass(frozen=True)
class TetrahedralComplex(StandardRegion):
    """Finite union of three-dimensional tetrahedra."""

    tetrahedra: tuple[Simplex, ...]

    def __init__(self, tetrahedra: Sequence[Simplex | Sequence[Sequence[object]]]):
        items = []
        for tetrahedron in tetrahedra:
            simplex = tetrahedron if isinstance(tetrahedron, Simplex) else Simplex(tetrahedron)
            if len(simplex.vertices) != 4 or simplex.ambient_dimension() != 3:
                raise ValueError("tetrahedra must have four three-dimensional vertices")
            items.append(simplex)
        object.__setattr__(self, "tetrahedra", tuple(items))

    def dimension(self) -> int:
        return max((tet.dimension() for tet in self.tetrahedra), default=-1)

    def ambient_dimension(self) -> int:
        return 3


@dataclass(frozen=True)
class Parallelogram(StandardRegion):
    origin: _PointData
    vectors: tuple[_PointData, _PointData]

    def __init__(self, origin: Sequence[object], vectors: Sequence[Sequence[object]]):
        if len(vectors) != 2:
            raise ValueError("a parallelogram requires two spanning vectors")
        origin_pt = _sympify_point(origin)
        vecs = tuple(_sympify_point(v) for v in vectors)
        if any(len(vec) != len(origin_pt) for vec in vecs):
            raise ValueError("parallelogram vectors must match the origin dimension")
        object.__setattr__(self, "origin", origin_pt)
        object.__setattr__(self, "vectors", vecs)

    def dimension(self) -> int:
        return _vector_rank(self.vectors)

    def ambient_dimension(self) -> int:
        return len(self.origin)


@dataclass(frozen=True)
class _ParallelepipedBase(StandardRegion):
    origin: _PointData
    vectors: tuple[_PointData, ...]

    def __init__(self, origin: Sequence[object], vectors: Sequence[Sequence[object]]):
        origin_pt = _sympify_point(origin)
        vecs = tuple(_sympify_point(v) for v in vectors)
        if any(len(vec) != len(origin_pt) for vec in vecs):
            raise ValueError("parallelepiped vectors must match the origin dimension")
        object.__setattr__(self, "origin", origin_pt)
        object.__setattr__(self, "vectors", vecs)

    def dimension(self) -> int:
        return _vector_rank(self.vectors)

    def ambient_dimension(self) -> int:
        return len(self.origin)


@dataclass(frozen=True)
class Parallelepiped(_ParallelepipedBase):
    """Canonical affine image of a unit box."""

    def __init__(
        self,
        origin: Sequence[object] | _ParallelepipedBase,
        vectors: Sequence[Sequence[object]] | None = None,
    ):
        if isinstance(origin, _ParallelepipedBase):
            if vectors is not None:
                raise TypeError("vectors must be omitted when converting a _ParallelepipedBase")
            vectors = origin.vectors
            origin = origin.origin
        if vectors is None:
            raise TypeError("vectors are required")
        super().__init__(origin, vectors)
        if self.dimension() != len(self.vectors):
            raise ValueError("parallelepiped spanning vectors must be linearly independent")

    def affine_hull(self):
        return AffineSpace(self.origin, self.vectors)

    @classmethod
    def from_region(cls, region: _ParallelepipedBase) -> Parallelepiped:
        return cls(region)


@dataclass(frozen=True)
class _PrismBase(StandardRegion):
    base: _PolygonBase | _SimplexBase
    vector: _PointData

    def __init__(
        self,
        base: _PolygonBase | _SimplexBase | Sequence[Sequence[object]],
        vector: Sequence[object],
    ):
        base_obj = base if isinstance(base, (_PolygonBase, _SimplexBase)) else _PolygonBase(base)
        vec = _sympify_point(vector)
        if len(vec) != base_obj.ambient_dimension():
            raise ValueError("prism vector must match the base ambient dimension")
        object.__setattr__(self, "base", base_obj)
        object.__setattr__(self, "vector", vec)

    def dimension(self) -> int:
        if isinstance(self.base, _PolygonBase):
            points = self.base.vertices
        else:
            points = self.base.vertices
        anchor = points[0]
        span = tuple(tuple(sp.Matrix(point) - sp.Matrix(anchor)) for point in points[1:]) + (
            self.vector,
        )
        return _vector_rank(span)

    def ambient_dimension(self) -> int:
        return len(self.vector)


@dataclass(frozen=True)
class _PyramidBase(StandardRegion):
    base: _PolygonBase | _SimplexBase
    apex: _PointData

    def __init__(
        self,
        base: _PolygonBase | _SimplexBase | Sequence[Sequence[object]],
        apex: Sequence[object],
    ):
        base_obj = base if isinstance(base, (_PolygonBase, _SimplexBase)) else _PolygonBase(base)
        apex_pt = _sympify_point(apex)
        if len(apex_pt) != base_obj.ambient_dimension():
            raise ValueError("pyramid apex must match the base ambient dimension")
        object.__setattr__(self, "base", base_obj)
        object.__setattr__(self, "apex", apex_pt)

    def dimension(self) -> int:
        points = self.base.vertices
        anchor = points[0]
        span = tuple(tuple(sp.Matrix(point) - sp.Matrix(anchor)) for point in points[1:]) + (
            tuple(sp.Matrix(self.apex) - sp.Matrix(anchor)),
        )
        return _vector_rank(span)

    def ambient_dimension(self) -> int:
        return len(self.apex)


def _positive_length(value: object, *, label: str) -> sp.Expr:
    expr = sp.sympify(value)
    if _provably_nonpositive(expr):
        raise ValueError(f"{label} must be positive")
    return expr


def _center3(center: Sequence[object]) -> _PointData:
    point = _sympify_point(center)
    if len(point) != 3:
        raise ValueError("solid center must be three-dimensional")
    return point


def _translate_scale_vertices(
    vertices: Sequence[Sequence[object]],
    center: Sequence[object],
    scale: sp.Expr,
) -> tuple[_PointData, ...]:
    c = _center3(center)
    return tuple(
        tuple(sp.simplify(c[i] + scale * sp.sympify(vertex[i])) for i in range(3))
        for vertex in vertices
    )


def Cube(
    center: Sequence[object] = (0, 0, 0),
    side: object = 1,
) -> Parallelepiped:
    """Return an axis-aligned cube as a canonical :class:`Parallelepiped`."""

    c = _center3(center)
    a = _positive_length(side, label="cube side length")
    half = a / 2
    origin = tuple(sp.simplify(x - half) for x in c)
    vectors = ((a, 0, 0), (0, a, 0), (0, 0, a))
    return Parallelepiped(origin, vectors)


def Tetrahedron(
    vertices: Sequence[Sequence[object]] | None = None,
    *,
    center: Sequence[object] = (0, 0, 0),
    edge: object = 1,
) -> Simplex:
    """Return a tetrahedron as a canonical :class:`Simplex`.

    With explicit ``vertices`` those vertices are used directly. Otherwise a
    regular tetrahedron centered at ``center`` with the requested edge length
    is constructed.
    """

    if vertices is not None:
        if center != (0, 0, 0) or edge != 1:
            raise TypeError("center and edge are unavailable with explicit vertices")
        if len(vertices) != 4:
            raise ValueError("a tetrahedron requires four vertices")
        simplex = Simplex(vertices)
        if simplex.ambient_dimension() != 3:
            raise ValueError("a tetrahedron requires three-dimensional vertices")
        return simplex
    a = _positive_length(edge, label="tetrahedron edge length")
    base = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    return Simplex(_translate_scale_vertices(base, center, a / (2 * sp.sqrt(2))))


def Octahedron(
    center: Sequence[object] = (0, 0, 0),
    edge: object = 1,
) -> Polytope:
    """Return a regular octahedron as a canonical :class:`Polytope`."""

    a = _positive_length(edge, label="octahedron edge length")
    base = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    return Polytope(_translate_scale_vertices(base, center, a / sp.sqrt(2)))


def Icosahedron(
    center: Sequence[object] = (0, 0, 0),
    edge: object = 1,
) -> Polytope:
    """Return a regular icosahedron as a canonical :class:`Polytope`."""

    a = _positive_length(edge, label="icosahedron edge length")
    phi = (1 + sp.sqrt(5)) / 2
    base = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            base.extend(((0, s1, s2 * phi), (s1, s2 * phi, 0), (s2 * phi, 0, s1)))
    return Polytope(_translate_scale_vertices(base, center, a / 2))


def Dodecahedron(
    center: Sequence[object] = (0, 0, 0),
    edge: object = 1,
) -> Polytope:
    """Return a regular dodecahedron as a canonical :class:`Polytope`."""

    a = _positive_length(edge, label="dodecahedron edge length")
    phi = (1 + sp.sqrt(5)) / 2
    inv = 1 / phi
    base = [(x, y, z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            base.extend(((0, s1 * inv, s2 * phi), (s1 * inv, s2 * phi, 0), (s2 * phi, 0, s1 * inv)))
    return Polytope(_translate_scale_vertices(base, center, a * phi / 2))


def Prism(
    base: Sequence[Sequence[object]] | StandardRegion,
    vector: Sequence[object],
) -> Polytope:
    """Extrude a vertex-defined base by ``vector`` and return a polytope."""

    if isinstance(base, (Polytope, _SimplexBase, _PolygonBase)):
        vertices = base.vertices
    else:
        vertices = tuple(_sympify_point(v) for v in base)
    if not vertices:
        raise ValueError("a prism requires a nonempty base")
    _validate_same_dimension(vertices, label="prism base vertices")
    vec = _sympify_point(vector)
    if len(vec) != len(vertices[0]):
        raise ValueError("prism vector must match the base ambient dimension")
    if all(_provably_zero(x) for x in vec):
        raise ValueError("prism vector must be nonzero")
    top = tuple(tuple(sp.simplify(x + d) for x, d in zip(v, vec, strict=True)) for v in vertices)
    return Polytope((*vertices, *top))


def Pyramid(
    base: Sequence[Sequence[object]] | StandardRegion,
    apex: Sequence[object],
) -> Polytope:
    """Join a vertex-defined base to an apex and return a polytope."""

    if isinstance(base, (Polytope, _SimplexBase, _PolygonBase)):
        vertices = base.vertices
    else:
        vertices = tuple(_sympify_point(v) for v in base)
    if not vertices:
        raise ValueError("a pyramid requires a nonempty base")
    _validate_same_dimension(vertices, label="pyramid base vertices")
    apex_pt = _sympify_point(apex)
    if len(apex_pt) != len(vertices[0]):
        raise ValueError("pyramid apex must match the base ambient dimension")
    result = Polytope((*vertices, apex_pt))
    if result.dimension() <= _affine_dimension(vertices):
        raise ValueError("pyramid apex must lie outside the base affine hull")
    return result


@dataclass(frozen=True)
class Hexahedron(Polytope):
    """Convex three-dimensional polytope with eight vertices and six quadrilateral facets."""

    def __init__(self, vertices: Sequence[Sequence[object]]):
        if len(vertices) != 8:
            raise ValueError("a hexahedron requires eight vertices")
        super().__init__(vertices)
        if self.ambient_dimension() != 3 or self.dimension() != 3:
            raise ValueError("a hexahedron requires full-dimensional three-dimensional vertices")
        if len(self.vertices) != 8:
            raise ValueError("a hexahedron requires eight distinct vertices")
        if not self.validate_topology(vertex_count=8, facet_count=6, facet_sizes=(4,) * 6):
            raise ValueError("vertices must form a convex hexahedron with six quadrilateral facets")

    @property
    def face_vertex_indices(self) -> tuple[tuple[int, ...], ...]:
        """Return the six quadrilateral facets as indices into :attr:`vertices`."""

        incidence = self.incidence()
        return tuple(tuple(facet.vertex_indices) for facet in incidence.facets)


@dataclass(frozen=True)
class Zonotope(_ParallelepipedBase):
    """Minkowski sum of line segments with possibly dependent generators."""
