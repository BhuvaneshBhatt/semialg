# ruff: noqa: F403, F405
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import sympy as sp

from ._standard_region_geometry import (
    _identity_directions,
    _PointData,
    _positive_definite_conditions,
    _provably_nonpositive,
    _provably_zero,
    _simplex_face_measure,
    _sphere_through_points,
    _sympify_point,
    _sympify_square_matrix,
    _validate_interval,
    _validate_nonnegative,
)
from ._standard_regions_base import *
from ._standard_regions_polyhedral import *


@dataclass(frozen=True)
class BallRegion(StandardRegion):
    center: _PointData
    radius: sp.Expr

    def __init__(self, center: Sequence[object], radius: object = 1, *, assumptions=None):
        radius_expr = sp.sympify(radius)
        _validate_nonnegative(radius_expr, label="radius")
        from .geometry_validity import validate_geometry_conditions

        validate_geometry_conditions((radius_expr >= 0,), assumptions, label="ball/sphere")
        object.__setattr__(self, "center", _sympify_point(center))
        object.__setattr__(self, "radius", radius_expr)

    @property
    def conditions(self) -> tuple[sp.Expr, ...]:
        condition = sp.simplify(self.radius >= 0)
        return () if condition is sp.true else (condition,)

    def dimension(self) -> int:
        return 0 if _provably_zero(self.radius) else len(self.center)

    def ambient_dimension(self) -> int:
        return len(self.center)


@dataclass(frozen=True)
class SphereRegion(BallRegion):
    def __init__(self, center: Sequence[object], radius: object = 1, *, assumptions=None):
        super().__init__(center, radius, assumptions=assumptions)

    def dimension(self) -> int:
        if _provably_zero(self.radius):
            return 0
        return max(0, len(self.center) - 1)


@dataclass(frozen=True)
class Ball(BallRegion):
    """Canonical closed ball in an arbitrary-dimensional Euclidean space."""

    def __init__(self, center: Sequence[object], radius: object = 1, *, assumptions=None):
        super().__init__(center, radius, assumptions=assumptions)

    def affine_hull(self):
        if _provably_zero(self.radius):
            return Point(self.center).affine_hull()
        return AffineSpace(self.center, _identity_directions(self.ambient_dimension()))

    @property
    def shape_matrix(self) -> sp.ImmutableMatrix:
        """Return the quadratic shape matrix ``radius**2 * I``."""

        return sp.ImmutableMatrix.eye(self.ambient_dimension()) * self.radius**2

    @classmethod
    def from_region(cls, region: BallRegion) -> Ball:
        """Convert an existing filled ball region to the canonical type."""

        if isinstance(region, SphereRegion):
            raise TypeError("a sphere boundary cannot be converted to Ball")
        return cls(region.center, region.radius)


@dataclass(frozen=True)
class Sphere(SphereRegion):
    """Canonical sphere boundary in an arbitrary-dimensional Euclidean space."""

    def __init__(self, center: Sequence[object], radius: object = 1, *, assumptions=None):
        super().__init__(center, radius, assumptions=assumptions)

    def affine_hull(self):
        if _provably_zero(self.radius):
            return Point(self.center).affine_hull()
        return AffineSpace(self.center, _identity_directions(self.ambient_dimension()))

    @property
    def shape_matrix(self) -> sp.ImmutableMatrix:
        """Return the quadratic shape matrix ``radius**2 * I``."""

        return sp.ImmutableMatrix.eye(self.ambient_dimension()) * self.radius**2

    @classmethod
    def from_region(cls, region: SphereRegion) -> Sphere:
        """Convert an existing sphere region to the canonical type."""

        return cls(region.center, region.radius)

    @classmethod
    def through(cls, *points: Sequence[object]) -> Sphere:
        """Construct the unique sphere through ``n + 1`` points in ``R^n``."""

        center, radius = _sphere_through_points(points)
        return cls(center, radius)

    @classmethod
    def circumscribed(cls, simplex: Simplex | SimplexRegion) -> Sphere:
        """Return the circumsphere of a full-dimensional simplex."""

        region = simplex if isinstance(simplex, Simplex) else Simplex(simplex)
        if region.dimension() != region.ambient_dimension():
            raise ValueError("circumsphere requires a full-dimensional simplex")
        return cls.through(*region.vertices)

    @classmethod
    def inscribed(cls, simplex: Simplex | SimplexRegion) -> Sphere:
        """Return the insphere of a full-dimensional simplex."""

        region = simplex if isinstance(simplex, Simplex) else Simplex(simplex)
        n = region.dimension()
        if n == 0 or n != region.ambient_dimension():
            raise ValueError("insphere requires a positive-dimensional full-dimensional simplex")
        vertices = region.vertices
        weights = tuple(
            _simplex_face_measure(vertices[:i] + vertices[i + 1 :]) for i in range(len(vertices))
        )
        total = sp.simplify(sum(weights))
        center = tuple(
            sp.simplify(sum(weights[i] * vertices[i][j] for i in range(len(vertices))) / total)
            for j in range(region.ambient_dimension())
        )
        volume = _simplex_face_measure(vertices)
        radius = sp.simplify(n * volume / total)
        return cls(center, radius)


class Circle:
    """Convenience constructor namespace for canonical two-dimensional spheres."""

    def __new__(cls, center: Sequence[object], radius: object = 1) -> Sphere:
        point = _sympify_point(center)
        if len(point) != 2:
            raise ValueError("Circle requires a two-dimensional center")
        return Sphere(point, radius)

    @staticmethod
    def through(*points: Sequence[object]) -> Sphere:
        """Construct the unique circle through three noncollinear planar points."""

        if len(points) != 3:
            raise ValueError("Circle.through requires exactly three points")
        if any(len(_sympify_point(point)) != 2 for point in points):
            raise ValueError("Circle.through requires two-dimensional points")
        return Sphere.through(*points)


@dataclass(frozen=True)
class Ellipsoid(StandardRegion):
    """Canonical filled ellipsoid ``(x-c).Q**-1.(x-c) <= 1``.

    ``shape_matrix`` is the positive-definite matrix ``Q``.  Symbolic
    positive-definiteness conditions that cannot be decided at construction
    time are retained in ``construction_conditions``.
    """

    center: _PointData
    shape_matrix: sp.ImmutableMatrix
    construction_conditions: tuple[sp.Expr, ...] = field(default=(), compare=False)

    def __init__(self, center: Sequence[object], shape_matrix: object, *, assumptions=None):
        point = _sympify_point(center)
        if not point:
            raise ValueError("an ellipsoid requires a nonempty center")
        matrix = _sympify_square_matrix(
            shape_matrix, dimension=len(point), label="ellipsoid shape matrix"
        )
        conditions = _positive_definite_conditions(matrix)
        from .geometry_validity import validate_geometry_conditions

        validate_geometry_conditions(conditions, assumptions, label="ellipsoid")
        matrix = sp.ImmutableMatrix(
            matrix.rows, matrix.cols, lambda i, j: sp.simplify((matrix[i, j] + matrix[j, i]) / 2)
        )
        object.__setattr__(self, "center", point)
        object.__setattr__(self, "shape_matrix", matrix)
        object.__setattr__(self, "construction_conditions", conditions)

    def dimension(self) -> int:
        return self.ambient_dimension()

    def ambient_dimension(self) -> int:
        return len(self.center)

    def affine_hull(self):
        return AffineSpace(self.center, _identity_directions(self.ambient_dimension()))

    @classmethod
    def from_radii(cls, center: Sequence[object], radii: Sequence[object]) -> Ellipsoid:
        """Construct an axis-aligned ellipsoid from strictly positive semiaxes."""

        point = _sympify_point(center)
        values = tuple(sp.sympify(radius) for radius in radii)
        if len(values) != len(point):
            raise ValueError("ellipsoid radii must match the center dimension")
        conditions = []
        for radius in values:
            if _provably_nonpositive(radius):
                raise ValueError("ellipsoid radii must be positive")
            condition = sp.simplify(radius > 0)
            if condition is not sp.true:
                conditions.append(condition)
        result = cls(point, sp.diag(*(radius**2 for radius in values)))
        merged = tuple(dict.fromkeys((*conditions, *result.construction_conditions)))
        object.__setattr__(result, "construction_conditions", merged)
        return result


@dataclass(frozen=True)
class EllipsoidBoundary(StandardRegion):
    """Canonical ellipsoid boundary ``(x-c).Q**-1.(x-c) = 1``."""

    center: _PointData
    shape_matrix: sp.ImmutableMatrix
    construction_conditions: tuple[sp.Expr, ...] = field(default=(), compare=False)

    def __init__(self, center: Sequence[object], shape_matrix: object, *, assumptions=None):
        point = _sympify_point(center)
        if not point:
            raise ValueError("an ellipsoid boundary requires a nonempty center")
        matrix = _sympify_square_matrix(
            shape_matrix, dimension=len(point), label="ellipsoid-boundary shape matrix"
        )
        conditions = _positive_definite_conditions(matrix)
        from .geometry_validity import validate_geometry_conditions

        validate_geometry_conditions(conditions, assumptions, label="ellipsoid boundary")
        matrix = sp.ImmutableMatrix(
            matrix.rows, matrix.cols, lambda i, j: sp.simplify((matrix[i, j] + matrix[j, i]) / 2)
        )
        object.__setattr__(self, "center", point)
        object.__setattr__(self, "shape_matrix", matrix)
        object.__setattr__(self, "construction_conditions", conditions)

    def dimension(self) -> int:
        return max(0, self.ambient_dimension() - 1)

    def ambient_dimension(self) -> int:
        return len(self.center)

    def affine_hull(self):
        return AffineSpace(self.center, _identity_directions(self.ambient_dimension()))


@dataclass(frozen=True)
class SphericalShellRegion(StandardRegion):
    center: _PointData
    inner_radius: sp.Expr
    outer_radius: sp.Expr

    def __init__(self, center: Sequence[object], radii: tuple[object, object]):
        inner = sp.sympify(radii[0])
        outer = sp.sympify(radii[1])
        _validate_nonnegative(inner, label="inner radius")
        _validate_nonnegative(outer, label="outer radius")
        _validate_interval(inner, outer, label="shell radii")
        object.__setattr__(self, "center", _sympify_point(center))
        object.__setattr__(self, "inner_radius", inner)
        object.__setattr__(self, "outer_radius", outer)

    @property
    def conditions(self) -> tuple[sp.Expr, ...]:
        conditions = (
            self.inner_radius >= 0,
            self.outer_radius >= 0,
            self.inner_radius <= self.outer_radius,
        )
        return tuple(cond for item in conditions if (cond := sp.simplify(item)) is not sp.true)

    def dimension(self) -> int:
        if _provably_zero(self.outer_radius - self.inner_radius):
            if _provably_zero(self.outer_radius):
                return 0
            return max(0, len(self.center) - 1)
        return len(self.center)

    def ambient_dimension(self) -> int:
        return len(self.center)


@dataclass(frozen=True)
class CylinderRegion(StandardRegion):
    start: _PointData
    end: _PointData
    radius: sp.Expr = sp.Integer(1)

    def __init__(self, start: Sequence[object], end: Sequence[object], radius: object = 1):
        start_pt = _sympify_point(start)
        end_pt = _sympify_point(end)
        if len(start_pt) != len(end_pt):
            raise ValueError("cylinder endpoints must have the same dimension")
        radius_expr = sp.sympify(radius)
        _validate_nonnegative(radius_expr, label="radius")
        object.__setattr__(self, "start", start_pt)
        object.__setattr__(self, "end", end_pt)
        object.__setattr__(self, "radius", radius_expr)

    @property
    def conditions(self) -> tuple[sp.Expr, ...]:
        radius_condition = sp.simplify(self.radius >= 0)
        inherited = tuple(getattr(self, "construction_conditions", ()))
        extra = () if radius_condition is sp.true else (radius_condition,)
        return tuple(dict.fromkeys((*extra, *inherited)))

    def dimension(self) -> int:
        axis_zero = all(_provably_zero(b - a) for a, b in zip(self.start, self.end, strict=True))
        if not _provably_zero(self.radius):
            return len(self.start)
        return 0 if axis_zero else 1

    def ambient_dimension(self) -> int:
        return len(self.start)


@dataclass(frozen=True)
class ConeRegion(CylinderRegion):
    def __init__(self, start: Sequence[object], end: Sequence[object], radius: object = 1):
        super().__init__(start, end, radius)


@dataclass(frozen=True)
class Cylinder(CylinderRegion):
    """Canonical flat-ended right circular cylinder.

    ``start`` and ``end`` are the centers of the two end caps.  The axis must
    be nondegenerate; a zero radius is allowed and represents the closed axis
    segment.
    """

    construction_conditions: tuple[sp.Expr, ...] = field(default=(), compare=False)

    def __init__(self, start: Sequence[object], end: Sequence[object], radius: object = 1):
        super().__init__(start, end, radius)
        if self.ambient_dimension() < 2:
            raise ValueError("Cylinder requires ambient dimension at least two")
        axis_sq = sp.simplify(sum((b - a) ** 2 for a, b in zip(self.start, self.end, strict=True)))
        if _provably_zero(axis_sq):
            raise ValueError("cylinder axis endpoints must be distinct")
        condition = sp.simplify(axis_sq > 0)
        object.__setattr__(
            self, "construction_conditions", () if condition is sp.true else (condition,)
        )

    @property
    def axis(self) -> _PointData:
        """Return the axis direction from ``start`` to ``end``."""

        return tuple(b - a for a, b in zip(self.start, self.end, strict=True))

    @property
    def height(self) -> sp.Expr:
        """Return the Euclidean length of the cylinder axis."""

        return sp.sqrt(sp.simplify(sum(value**2 for value in self.axis)))

    def affine_hull(self):
        if _provably_zero(self.radius):
            return Line(self.start, self.axis)
        return AffineSpace(self.start, _identity_directions(self.ambient_dimension()))

    @classmethod
    def from_region(cls, region: CylinderRegion) -> Cylinder:
        """Convert an existing nondegenerate cylinder region."""

        if isinstance(region, ConeRegion):
            raise TypeError("cannot convert ConeRegion to Cylinder")
        return cls(region.start, region.end, region.radius)


@dataclass(frozen=True)
class Cone(ConeRegion):
    """Canonical right circular cone with base at ``start`` and apex at ``end``."""

    construction_conditions: tuple[sp.Expr, ...] = field(default=(), compare=False)

    def __init__(self, start: Sequence[object], end: Sequence[object], radius: object = 1):
        super().__init__(start, end, radius)
        if self.ambient_dimension() < 2:
            raise ValueError("Cone requires ambient dimension at least two")
        axis_sq = sp.simplify(sum((b - a) ** 2 for a, b in zip(self.start, self.end, strict=True)))
        if _provably_zero(axis_sq):
            raise ValueError("cone axis endpoints must be distinct")
        condition = sp.simplify(axis_sq > 0)
        object.__setattr__(
            self, "construction_conditions", () if condition is sp.true else (condition,)
        )

    @property
    def axis(self) -> _PointData:
        """Return the axis direction from the base center to the apex."""

        return tuple(b - a for a, b in zip(self.start, self.end, strict=True))

    @property
    def height(self) -> sp.Expr:
        """Return the Euclidean base-to-apex height."""

        return sp.sqrt(sp.simplify(sum(value**2 for value in self.axis)))

    def affine_hull(self):
        if _provably_zero(self.radius):
            return Line(self.start, self.axis)
        return AffineSpace(self.start, _identity_directions(self.ambient_dimension()))

    @classmethod
    def from_region(cls, region: ConeRegion) -> Cone:
        """Convert an existing nondegenerate cone region."""

        return cls(region.start, region.end, region.radius)


@dataclass(frozen=True)
class Torus(StandardRegion):
    """Canonical ring or horn torus surface in three-dimensional space."""

    center: _PointData
    major_radius: sp.Expr
    minor_radius: sp.Expr
    construction_conditions: tuple[sp.Expr, ...] = field(default=(), compare=False)

    def __init__(self, center: Sequence[object], major_radius: object, minor_radius: object):
        point = _sympify_point(center)
        if len(point) != 3:
            raise ValueError("Torus requires a three-dimensional center")
        major = sp.sympify(major_radius)
        minor = sp.sympify(minor_radius)
        if _provably_nonpositive(major) or _provably_nonpositive(minor):
            raise ValueError("torus radii must be positive")
        gap = sp.simplify(major - minor)
        if gap.is_negative is True:
            raise ValueError("major radius must be at least the minor radius")
        conditions = []
        for condition in (major > 0, minor > 0, major >= minor):
            cond = sp.simplify(condition)
            if cond is not sp.true:
                conditions.append(cond)
        object.__setattr__(self, "center", point)
        object.__setattr__(self, "major_radius", major)
        object.__setattr__(self, "minor_radius", minor)
        object.__setattr__(self, "construction_conditions", tuple(dict.fromkeys(conditions)))

    def dimension(self) -> int:
        return 2

    def ambient_dimension(self) -> int:
        return 3

    def affine_hull(self):
        return AffineSpace(self.center, _identity_directions(3))

    @property
    def inner_radius(self) -> sp.Expr:
        """Return the distance from the axis to the inner equator."""
        return sp.simplify(self.major_radius - self.minor_radius)

    @property
    def outer_radius(self) -> sp.Expr:
        """Return the distance from the axis to the outer equator."""
        return sp.simplify(self.major_radius + self.minor_radius)


@dataclass(frozen=True)
class FilledTorus(Torus):
    """Canonical solid ring or horn torus in three-dimensional space."""

    def dimension(self) -> int:
        return 3


@dataclass(frozen=True)
class StadiumRegion(StandardRegion):
    start: _PointData
    end: _PointData
    radius: sp.Expr = sp.Integer(1)

    def __init__(self, start: Sequence[object], end: Sequence[object], radius: object = 1):
        start_pt = _sympify_point(start)
        end_pt = _sympify_point(end)
        if len(start_pt) != len(end_pt):
            raise ValueError("region endpoints must have the same dimension")
        if len(start_pt) != 2:
            raise ValueError("StadiumRegion requires two-dimensional endpoints")
        radius_expr = sp.sympify(radius)
        _validate_nonnegative(radius_expr, label="radius")
        object.__setattr__(self, "start", start_pt)
        object.__setattr__(self, "end", end_pt)
        object.__setattr__(self, "radius", radius_expr)

    @property
    def conditions(self) -> tuple[sp.Expr, ...]:
        condition = sp.simplify(self.radius >= 0)
        return () if condition is sp.true else (condition,)

    def dimension(self) -> int:
        axis_zero = all(_provably_zero(b - a) for a, b in zip(self.start, self.end, strict=True))
        if not _provably_zero(self.radius):
            return 2
        return 0 if axis_zero else 1

    def ambient_dimension(self) -> int:
        return 2


@dataclass(frozen=True)
class CapsuleRegion(StadiumRegion):
    def __init__(self, start: Sequence[object], end: Sequence[object], radius: object = 1):
        start_pt = _sympify_point(start)
        end_pt = _sympify_point(end)
        if len(start_pt) != len(end_pt):
            raise ValueError("region endpoints must have the same dimension")
        radius_expr = sp.sympify(radius)
        _validate_nonnegative(radius_expr, label="radius")
        object.__setattr__(self, "start", start_pt)
        object.__setattr__(self, "end", end_pt)
        object.__setattr__(self, "radius", radius_expr)

    def dimension(self) -> int:
        axis_zero = all(_provably_zero(b - a) for a, b in zip(self.start, self.end, strict=True))
        if not _provably_zero(self.radius):
            return len(self.start)
        return 0 if axis_zero else 1

    def ambient_dimension(self) -> int:
        return len(self.start)
