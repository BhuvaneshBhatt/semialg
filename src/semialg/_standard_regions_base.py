from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ._standard_region_geometry import (
    _independent_vectors,
    _PointData,
    _provably_zero,
    _sympify_point,
    _validate_same_dimension,
    _vector_rank,
)


class Geometry:
    """Common protocol for explicit geometric objects.

    Geometry objects preserve structural information while providing a single
    lowering path to semialgebraic formulas.  Shape-specific capabilities may
    be exposed by subclasses without becoming requirements of the base
    protocol.
    """

    def dimension(self) -> int:
        """Return the intrinsic dimension of the geometry."""

        raise NotImplementedError

    def ambient_dimension(self) -> int:
        """Return the dimension of the containing coordinate space."""

        raise NotImplementedError

    @property
    def intrinsic_dimension(self) -> int:
        """Return the intrinsic dimension using property-style access."""

        return self.dimension()

    def as_semialgebraic_region(self, variables=None):
        """Return the unified symbolic-region representation of this geometry."""

        from .symbolic_regions import as_semialgebraic_region

        return as_semialgebraic_region(self, variables)

    def as_formula(self, variables=None, *, eliminate: bool = False) -> sp.Expr:
        """Return a semialgebraic membership formula for this geometry.

        By default the formula preserves any existential parameters introduced
        by structural lowering.  Set ``eliminate=True`` to request the
        quantifier-free formula used by CAD and decision procedures.
        """

        region = self.as_semialgebraic_region(variables)
        return region.quantifier_free_formula() if eliminate else region.formula

    def contains(self, point, *, strategy: str = "auto") -> bool:
        """Return whether an exact point belongs to the geometry."""

        return self.as_semialgebraic_region().contains(point, strategy=strategy)

    def boundary(self, variables=None, *, strategy: str | None = None):
        """Return the Euclidean boundary as a :class:`SemialgebraicRegion`."""

        from .regions.operations import region_boundary

        return region_boundary(self, variables, strategy=strategy)

    def subset_of(self, other, *, strategy: str | None = None) -> bool:
        """Return whether this geometry is contained in ``other``."""

        from .reasoning_regions import region_subset

        return region_subset(self, other, strategy=strategy)

    def disjoint_from(self, other, *, strategy: str | None = None) -> bool:
        """Return whether this geometry is disjoint from ``other``."""

        from .reasoning_regions import region_disjoint

        return region_disjoint(self, other, strategy=strategy)

    def equals_region(self, other, *, strategy: str | None = None) -> bool:
        """Return whether this geometry and ``other`` define the same set."""

        from .reasoning_regions import region_equal

        return region_equal(self, other, strategy=strategy)

    def intersection(self, other):
        """Return the exact intersection with another semialgebraic region."""

        from .regions.operations import region_intersection

        return region_intersection(self, other)

    def product(self, other, *, variables=None):
        """Return the Cartesian product with another semialgebraic region."""

        from .regions.operations import region_product

        return region_product(self, other, variables=variables)

    def minkowski_sum(self, other, variables=None):
        """Return the exact Minkowski sum with another semialgebraic region."""

        from .derived_geometry import minkowski_sum

        return minkowski_sum(self, other, variables)

    def image(self, mapping, variables=None):
        """Return the exact image under an affine or polynomial mapping."""

        from .region_transformations import region_image

        return region_image(self, mapping, variables)

    def preimage(self, mapping, variables=None, *, target_variables=None):
        """Return the exact preimage under an affine or polynomial mapping."""

        from .region_transformations import region_preimage

        return region_preimage(self, mapping, variables, target_variables=target_variables)

    def affine_hull(self):
        """Return a structural affine hull when a subclass provides one."""

        raise NotImplementedError(f"{type(self).__name__} does not expose a structural affine hull")

    @property
    def conditions(self) -> tuple[sp.Expr, ...]:
        """Return symbolic conditions required for this geometry to be valid.

        Constructors reject conditions that are already provably false, while
        undecidable symbolic requirements remain explicit here.
        """

        return tuple(getattr(self, "construction_conditions", ()))

    def is_valid(self, assumptions=None, *, strategy: str | None = None) -> bool | None:
        """Decide validity under assumptions when possible.

        ``True`` means all :attr:`conditions` are proved, ``False`` means at
        least one is disproved, and ``None`` means the available symbolic
        information is insufficient.
        """

        from .geometry_validity import validity_status

        return validity_status(self.conditions, assumptions, strategy=strategy)

    def sample_point(self, variables=None, *, strategy: str = "auto", **kwargs):
        """Return one certified representative point as an ambient coordinate tuple."""

        from .sampling import sample_point

        symbolic = self.as_semialgebraic_region(variables)
        point = sample_point(
            symbolic.quantifier_free_formula(), symbolic.variables, strategy=strategy, **kwargs
        )
        if point is None:
            return None
        return tuple(point[var] for var in symbolic.variables)

    def sample_points(self, variables=None, *, count: int = 1, strategy: str = "auto", **kwargs):
        """Return certified representative points as ambient coordinate tuples."""

        from .sampling import sample_points

        symbolic = self.as_semialgebraic_region(variables)
        points = sample_points(
            symbolic.quantifier_free_formula(),
            symbolic.variables,
            count=count,
            strategy=strategy,
            **kwargs,
        )
        return tuple(tuple(point[var] for var in symbolic.variables) for point in points)

    def random_point(self, *, seed=None, distribution: str = "uniform", precision: int = 30):
        """Return a random point, uniform in intrinsic measure when supported."""

        from .region_sampling import random_point

        return random_point(self, seed=seed, distribution=distribution, precision=precision)

    def random_points(
        self, count: int, *, seed=None, distribution: str = "uniform", precision: int = 30
    ):
        """Return random points, uniform in intrinsic measure when supported."""

        from .region_sampling import random_points

        return random_points(
            self, count=count, seed=seed, distribution=distribution, precision=precision
        )

    def measure(self, *, measure_dimension: object = "intrinsic", return_result: bool = False):
        """Return exact ambient or intrinsic measure of this geometry."""

        from .region_measurement import region_measure

        return region_measure(
            self, measure_dimension=measure_dimension, return_result=return_result
        )

    def centroid(self, *, measure_dimension: object = "intrinsic", return_result: bool = False):
        """Return the centroid with respect to the requested uniform measure."""

        from .region_measurement import geometry_centroid

        return geometry_centroid(
            self, measure_dimension=measure_dimension, return_result=return_result
        )

    def transform(self, matrix, offset=None):
        """Return the exact affine image ``A*x + b`` of this geometry.

        Canonical structure is preserved whenever the transformed set has a
        supported canonical representation; otherwise an exact
        :class:`TransformedRegion` is returned.
        """

        from .region_transformations import affine_image

        return affine_image(self, matrix, offset)


class StandardRegion(Geometry):
    """Base class for explicit standard geometry supported by semialg."""

    def bounded_parametric_cover(self, variables=None, bounds=None):
        """Return a certified bounded parametric cover when structurally available."""

        from .parametric_geometry import bounded_parametric_cover

        return bounded_parametric_cover(self, variables, bounds)

    def intrinsic_parametric_cover(self, variables=None, *, dimension=None, require_verified=True):
        """Return exact native charts for the region's intrinsic geometry."""

        from .parametric_geometry import intrinsic_parametric_cover

        return intrinsic_parametric_cover(
            self, variables, dimension=dimension, require_verified=require_verified
        )


@dataclass(frozen=True)
class PointRegion(StandardRegion):
    points: tuple[_PointData, ...]

    def __init__(self, points: Sequence[Sequence[object]] | Sequence[object]):
        if points and not isinstance(points[0], (list, tuple)):  # type: ignore[index]
            pts = (_sympify_point(points),)  # type: ignore[arg-type]
        else:
            pts = tuple(_sympify_point(p) for p in points)  # type: ignore[arg-type]
        _validate_same_dimension(pts, label="points")
        object.__setattr__(self, "points", pts)

    def dimension(self) -> int:
        return 0 if self.points else -1

    def ambient_dimension(self) -> int:
        return len(self.points[0]) if self.points else 0


@dataclass(frozen=True)
class Point(PointRegion):
    """A single point in an affine coordinate space."""

    def __init__(self, coordinates: Sequence[object]):
        super().__init__(coordinates)

    @property
    def coordinates(self) -> _PointData:
        """Return the point coordinates."""

        return self.points[0]

    def affine_hull(self):
        """Return the zero-dimensional affine space containing the point."""

        return AffineSpace(self.coordinates, ())


def _validate_vector(vector: Sequence[object], *, label: str) -> _PointData:
    values = _sympify_point(vector)
    if not values:
        raise ValueError(f"{label} must not be empty")
    if all(_provably_zero(value) for value in values):
        raise ValueError(f"{label} must be nonzero")
    return values


@dataclass(frozen=True)
class AffineSpace(StandardRegion):
    """Affine space represented by one point and independent directions."""

    point: _PointData
    directions: tuple[_PointData, ...]

    def __init__(self, point: Sequence[object], directions: Sequence[Sequence[object]] = ()):
        anchor = _sympify_point(point)
        dirs = tuple(_sympify_point(v) for v in directions)
        if not anchor:
            raise ValueError("an affine space requires a nonempty point")
        _validate_same_dimension((anchor, *dirs), label="affine-space data")
        basis = _independent_vectors(
            tuple(v for v in dirs if any(not _provably_zero(x) for x in v))
        )
        object.__setattr__(self, "point", anchor)
        object.__setattr__(self, "directions", basis)

    def dimension(self) -> int:
        return len(self.directions)

    def ambient_dimension(self) -> int:
        return len(self.point)

    def affine_hull(self):
        return self


@dataclass(frozen=True)
class Hyperplane(StandardRegion):
    """Codimension-one affine space ``normal · (x - point) = 0``."""

    normal: _PointData
    point: _PointData

    def __init__(self, normal: Sequence[object], point: Sequence[object]):
        n = _validate_vector(normal, label="hyperplane normal")
        p = _sympify_point(point)
        _validate_same_dimension((n, p), label="hyperplane data")
        object.__setattr__(self, "normal", n)
        object.__setattr__(self, "point", p)

    def dimension(self) -> int:
        return self.ambient_dimension() - 1

    def ambient_dimension(self) -> int:
        return len(self.point)

    def affine_hull(self):
        matrix = sp.Matrix([self.normal])
        directions = tuple(tuple(v) for v in matrix.nullspace())
        return AffineSpace(self.point, directions)


@dataclass(frozen=True)
class HalfSpace(StandardRegion):
    """Closed half-space ``normal · (x - point) <= 0``."""

    normal: _PointData
    point: _PointData

    def __init__(self, normal: Sequence[object], point: Sequence[object]):
        n = _validate_vector(normal, label="half-space normal")
        p = _sympify_point(point)
        _validate_same_dimension((n, p), label="half-space data")
        object.__setattr__(self, "normal", n)
        object.__setattr__(self, "point", p)

    def dimension(self) -> int:
        return self.ambient_dimension()

    def ambient_dimension(self) -> int:
        return len(self.point)

    def affine_hull(self):
        return AffineSpace(self.point, sp.eye(self.ambient_dimension()).tolist())


@dataclass(frozen=True)
class Line(AffineSpace):
    """Infinite affine line through ``point`` in ``direction``."""

    def __init__(self, point: Sequence[object], direction: Sequence[object]):
        direction_data = _validate_vector(direction, label="line direction")
        super().__init__(point, (direction_data,))

    @property
    def direction(self) -> _PointData:
        return self.directions[0]


@dataclass(frozen=True)
class Ray(StandardRegion):
    """Closed ray ``point + t*direction`` for ``t >= 0``."""

    point: _PointData
    direction: _PointData

    def __init__(self, point: Sequence[object], direction: Sequence[object]):
        p = _sympify_point(point)
        d = _validate_vector(direction, label="ray direction")
        _validate_same_dimension((p, d), label="ray data")
        object.__setattr__(self, "point", p)
        object.__setattr__(self, "direction", d)

    def dimension(self) -> int:
        return 1

    def ambient_dimension(self) -> int:
        return len(self.point)

    def affine_hull(self):
        return Line(self.point, self.direction)


@dataclass(frozen=True)
class AffineHalfSpace(StandardRegion):
    """Half-space intrinsic to an affine subspace.

    ``directions`` span the boundary directions and ``inward`` selects the
    nonnegative half-direction from ``point``.
    """

    point: _PointData
    directions: tuple[_PointData, ...]
    inward: _PointData

    def __init__(
        self,
        point: Sequence[object],
        directions: Sequence[Sequence[object]],
        inward: Sequence[object],
    ):
        p = _sympify_point(point)
        dirs = tuple(_sympify_point(v) for v in directions)
        w = _validate_vector(inward, label="affine half-space inward direction")
        _validate_same_dimension((p, *dirs, w), label="affine half-space data")
        boundary_basis = _independent_vectors(
            tuple(v for v in dirs if any(not _provably_zero(x) for x in v))
        )
        if _vector_rank((*boundary_basis, w)) == len(boundary_basis):
            raise ValueError("inward direction must not lie in the boundary span")
        object.__setattr__(self, "point", p)
        object.__setattr__(self, "directions", boundary_basis)
        object.__setattr__(self, "inward", w)

    def dimension(self) -> int:
        return len(self.directions) + 1

    def ambient_dimension(self) -> int:
        return len(self.point)

    def affine_hull(self):
        return AffineSpace(self.point, (*self.directions, self.inward))
