from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .exact_arithmetic import compare_exact_reals
from .normalization import normalize_symbol_sequence

Point = tuple[sp.Expr, ...]


def _sympify_point(point: Sequence[object]) -> Point:
    return tuple(sp.sympify(v) for v in point)


def _validate_same_dimension(points: Sequence[Point], *, label: str) -> int:
    """Validate that all points/vectors use one ambient dimension."""

    if not points:
        return 0
    dim = len(points[0])
    if any(len(point) != dim for point in points[1:]):
        raise ValueError(f"{label} must all have the same dimension")
    return dim


def _validate_nonnegative(value: sp.Expr, *, label: str) -> None:
    """Reject values that are provably negative, allowing symbolic unknowns."""

    try:
        if compare_exact_reals(value, sp.Integer(0)) < 0:
            raise ValueError(f"{label} must be nonnegative") from None
    except (TypeError, ValueError, NotImplementedError):
        if sp.simplify(value).is_negative is True:
            raise ValueError(f"{label} must be nonnegative") from None


def _validate_interval(lower: sp.Expr, upper: sp.Expr, *, label: str = "bounds") -> None:
    """Reject an interval whose endpoint order is provably reversed."""

    try:
        reversed_order = compare_exact_reals(lower, upper) > 0
    except (TypeError, ValueError, NotImplementedError):
        reversed_order = sp.simplify(lower - upper).is_positive is True
    if reversed_order:
        raise ValueError(f"{label} have lower endpoint greater than upper endpoint")


def _provably_zero(value: sp.Expr) -> bool:
    """Return whether an exact symbolic quantity can be certified as zero."""

    simplified = sp.simplify(value)
    if simplified == 0 or simplified.is_zero is True:
        return True
    try:
        return compare_exact_reals(simplified, sp.Integer(0)) == 0
    except (TypeError, ValueError, NotImplementedError):
        return False


def _vector_rank(vectors: Sequence[Point]) -> int:
    """Return the exact/generic rank of a finite family of coordinate vectors."""

    if not vectors:
        return 0
    return int(sp.Matrix.hstack(*(sp.Matrix(vector) for vector in vectors)).rank())


def _affine_dimension(points: Sequence[Point]) -> int:
    """Return the affine dimension of a finite point set."""

    if len(points) < 2:
        return 0
    anchor = sp.Matrix(points[0])
    vectors = tuple(tuple(sp.Matrix(point) - anchor) for point in points[1:])
    return _vector_rank(vectors)


def _independent_vectors(vectors: Sequence[Point]) -> tuple[Point, ...]:
    """Return a deterministic maximal independent subfamily."""

    basis: list[Point] = []
    rank = 0
    for vector in vectors:
        candidate = (*basis, vector)
        new_rank = _vector_rank(candidate)
        if new_rank > rank:
            basis.append(vector)
            rank = new_rank
    return tuple(basis)


def _effective_parametric_data(
    parameters: Sequence[sp.Symbol],
    limits: Sequence[tuple[sp.Symbol, sp.Expr, sp.Expr]],
    mapping: Sequence[sp.Expr],
    assumptions: sp.Expr,
) -> tuple[
    tuple[sp.Symbol, ...],
    tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...],
    tuple[sp.Expr, ...],
    sp.Expr,
]:
    """Remove parameters whose integration interval is identically a point.

    Fixed values are substituted into all remaining bounds, the mapping, and
    assumptions.  Repeating to a fixed point handles chains such as ``v=u``
    followed by a fixed ``u`` limit.
    """

    fixed: dict[sp.Symbol, sp.Expr] = {}
    remaining = list(limits)
    changed = True
    while changed:
        changed = False
        kept: list[tuple[sp.Symbol, sp.Expr, sp.Expr]] = []
        for param, lower, upper in remaining:
            lo = sp.simplify(lower.subs(fixed))
            hi = sp.simplify(upper.subs(fixed))
            if _provably_zero(hi - lo):
                fixed[param] = lo
                fixed = {key: sp.simplify(value.subs(fixed)) for key, value in fixed.items()}
                changed = True
            else:
                kept.append((param, lo, hi))
        remaining = kept

    free = tuple(param for param in parameters if param not in fixed)
    final_limits = tuple(
        (param, sp.simplify(lower.subs(fixed)), sp.simplify(upper.subs(fixed)))
        for param, lower, upper in remaining
    )
    final_mapping = tuple(sp.simplify(expr.subs(fixed)) for expr in mapping)
    final_assumptions = sp.simplify(assumptions.subs(fixed))
    return free, final_limits, final_mapping, final_assumptions


def _orientation(a: Point, b: Point, c: Point) -> int:
    """Return the certified orientation sign of three exact planar points."""

    cross = sp.expand((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
    try:
        return compare_exact_reals(cross, sp.Integer(0))
    except (TypeError, ValueError, NotImplementedError) as exc:
        raise ValueError("polygon orientation could not be certified exactly") from exc


def _on_segment(a: Point, b: Point, p: Point) -> bool:
    if _orientation(a, b, p) != 0:
        return False
    for u, v, w in zip(a, b, p, strict=True):
        if compare_exact_reals(u, v) <= 0:
            lower, upper = u, v
        else:
            lower, upper = v, u
        if compare_exact_reals(lower, w) > 0 or compare_exact_reals(w, upper) > 0:
            return False
    return True


def _segments_intersect(a: Point, b: Point, c: Point, d: Point) -> bool:
    ab_c = _orientation(a, b, c)
    ab_d = _orientation(a, b, d)
    cd_a = _orientation(c, d, a)
    cd_b = _orientation(c, d, b)
    if ab_c * ab_d < 0 and cd_a * cd_b < 0:
        return True
    return (
        (ab_c == 0 and _on_segment(a, b, c))
        or (ab_d == 0 and _on_segment(a, b, d))
        or (cd_a == 0 and _on_segment(c, d, a))
        or (cd_b == 0 and _on_segment(c, d, b))
    )


def _validate_simple_polygon(vertices: Sequence[Point]) -> None:
    """Reject repeated vertices, zero edges, and self intersections."""

    if len(set(vertices)) != len(vertices):
        raise ValueError("polygon vertices must be distinct")
    count = len(vertices)
    for index in range(count):
        a, b = vertices[index], vertices[(index + 1) % count]
        if a == b:
            raise ValueError("polygon edges must have nonzero length")
        for other in range(index + 1, count):
            if other in {index, (index + 1) % count}:
                continue
            if index in {other, (other + 1) % count}:
                continue
            c, d = vertices[other], vertices[(other + 1) % count]
            if _segments_intersect(a, b, c, d):
                raise ValueError("polygon edges must not self-intersect")


def _signed_polygon_orientation(vertices: Sequence[Point]) -> int:
    twice_area = sp.expand(
        sum(
            a[0] * b[1] - b[0] * a[1]
            for a, b in zip(vertices, (*vertices[1:], vertices[0]), strict=True)
        )
    )
    try:
        orientation = compare_exact_reals(twice_area, sp.Integer(0))
    except (TypeError, ValueError, NotImplementedError) as exc:
        raise ValueError("polygon area could not be certified exactly") from exc
    if orientation == 0:
        raise ValueError("polygon vertices must enclose nonzero area")
    return orientation


def _point_in_triangle(point: Point, triangle: tuple[Point, Point, Point], sign: int) -> bool:
    a, b, c = triangle
    return all(
        sign * orient >= 0
        for orient in (
            _orientation(a, b, point),
            _orientation(b, c, point),
            _orientation(c, a, point),
        )
    )


def _triangulate_polygon(vertices: Sequence[Point]) -> tuple[tuple[Point, Point, Point], ...]:
    """Triangulate a certified simple polygon by exact ear clipping."""

    sign = _signed_polygon_orientation(vertices)
    remaining = list(range(len(vertices)))
    triangles: list[tuple[Point, Point, Point]] = []
    while len(remaining) > 3:
        ear_found = False
        for pos, current in enumerate(remaining):
            previous = remaining[pos - 1]
            following = remaining[(pos + 1) % len(remaining)]
            triangle = (vertices[previous], vertices[current], vertices[following])
            if sign * _orientation(*triangle) <= 0:
                continue
            if any(
                _point_in_triangle(vertices[index], triangle, sign)
                for index in remaining
                if index not in {previous, current, following}
            ):
                continue
            triangles.append(triangle)
            del remaining[pos]
            ear_found = True
            break
        if not ear_found:
            raise ValueError("polygon could not be triangulated as a simple polygon")
    triangles.append(tuple(vertices[index] for index in remaining))  # type: ignore[arg-type]
    return tuple(triangles)


class StandardRegion:
    """Base class for explicit region objects supported by semialg."""

    def dimension(self) -> int:
        raise NotImplementedError

    def ambient_dimension(self) -> int:
        raise NotImplementedError

    def as_semialgebraic_region(self, variables=None):
        """Return the unified symbolic-region representation of this region."""

        from .symbolic_regions import as_semialgebraic_region

        return as_semialgebraic_region(self, variables)

    def bounded_parametric_cover(self, variables=None, bounds=None):
        """Return a certified bounded parametric cover when structurally available."""

        from .parametric_geometry import bounded_parametric_cover

        return bounded_parametric_cover(self, variables, bounds)


@dataclass(frozen=True)
class PointRegion(StandardRegion):
    points: tuple[Point, ...]

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
class IntervalRegion(StandardRegion):
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
class BoxRegion(StandardRegion):
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
class SimplexRegion(StandardRegion):
    vertices: tuple[Point, ...]

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
class PolygonRegion(StandardRegion):
    vertices: tuple[Point, ...]

    def __init__(self, vertices: Sequence[Sequence[object]]):
        verts = tuple(_sympify_point(v) for v in vertices)
        if len(verts) > 1 and verts[0] == verts[-1]:
            verts = verts[:-1]
        if len(verts) < 3:
            raise ValueError("a polygon requires at least three vertices")
        if len({len(v) for v in verts}) != 1 or len(verts[0]) != 2:
            raise ValueError("PolygonRegion represents 2D polygons")
        _validate_simple_polygon(verts)
        _signed_polygon_orientation(verts)
        object.__setattr__(self, "vertices", verts)

    def dimension(self) -> int:
        return 2

    def ambient_dimension(self) -> int:
        return 2

    def triangulation(self) -> tuple[SimplexRegion, ...]:
        return tuple(SimplexRegion(triangle) for triangle in _triangulate_polygon(self.vertices))


@dataclass(frozen=True)
class TetrahedronRegion(SimplexRegion):
    def __init__(self, vertices: Sequence[Sequence[object]]):
        if len(vertices) != 4:
            raise ValueError("a tetrahedron requires four vertices")
        super().__init__(vertices)
        if self.ambient_dimension() != 3:
            raise ValueError("a tetrahedron requires three-dimensional vertices")


@dataclass(frozen=True)
class PolyhedronRegion(StandardRegion):
    tetrahedra: tuple[TetrahedronRegion, ...]

    def __init__(self, tetrahedra: Sequence[TetrahedronRegion | Sequence[Sequence[object]]]):
        tets = tuple(
            t if isinstance(t, TetrahedronRegion) else TetrahedronRegion(t) for t in tetrahedra
        )
        if any(tet.ambient_dimension() != 3 for tet in tets):
            raise ValueError("PolyhedronRegion tetrahedra must be three-dimensional")
        object.__setattr__(self, "tetrahedra", tets)

    def dimension(self) -> int:
        return max((tet.dimension() for tet in self.tetrahedra), default=-1)

    def ambient_dimension(self) -> int:
        return 3


@dataclass(frozen=True)
class ParallelogramRegion(StandardRegion):
    origin: Point
    vectors: tuple[Point, Point]

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
class ParallelepipedRegion(StandardRegion):
    origin: Point
    vectors: tuple[Point, ...]

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
class PrismRegion(StandardRegion):
    base: PolygonRegion | SimplexRegion
    vector: Point

    def __init__(
        self,
        base: PolygonRegion | SimplexRegion | Sequence[Sequence[object]],
        vector: Sequence[object],
    ):
        base_obj = base if isinstance(base, (PolygonRegion, SimplexRegion)) else PolygonRegion(base)
        vec = _sympify_point(vector)
        if len(vec) != base_obj.ambient_dimension():
            raise ValueError("prism vector must match the base ambient dimension")
        object.__setattr__(self, "base", base_obj)
        object.__setattr__(self, "vector", vec)

    def dimension(self) -> int:
        if isinstance(self.base, PolygonRegion):
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
class PyramidRegion(StandardRegion):
    base: PolygonRegion | SimplexRegion
    apex: Point

    def __init__(
        self,
        base: PolygonRegion | SimplexRegion | Sequence[Sequence[object]],
        apex: Sequence[object],
    ):
        base_obj = base if isinstance(base, (PolygonRegion, SimplexRegion)) else PolygonRegion(base)
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


@dataclass(frozen=True)
class BallRegion(StandardRegion):
    center: Point
    radius: sp.Expr

    def __init__(self, center: Sequence[object], radius: object = 1):
        radius_expr = sp.sympify(radius)
        _validate_nonnegative(radius_expr, label="radius")
        object.__setattr__(self, "center", _sympify_point(center))
        object.__setattr__(self, "radius", radius_expr)

    def dimension(self) -> int:
        return 0 if _provably_zero(self.radius) else len(self.center)

    def ambient_dimension(self) -> int:
        return len(self.center)


@dataclass(frozen=True)
class SphereRegion(BallRegion):
    def __init__(self, center: Sequence[object], radius: object = 1):
        super().__init__(center, radius)

    def dimension(self) -> int:
        if _provably_zero(self.radius):
            return 0
        return max(0, len(self.center) - 1)


@dataclass(frozen=True)
class SphericalShellRegion(StandardRegion):
    center: Point
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
    start: Point
    end: Point
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
class StadiumRegion(StandardRegion):
    start: Point
    end: Point
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


@dataclass(frozen=True)
class ParametricRegion(StandardRegion):
    parameters: tuple[sp.Symbol, ...]
    limits: tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]
    mapping: tuple[sp.Expr, ...]
    multiplicity: sp.Expr = sp.Integer(1)
    assumptions: sp.Expr = sp.true

    def __init__(
        self,
        parameters: Sequence[sp.Symbol | str],
        limits: Sequence[tuple[sp.Symbol | str, object, object]],
        mapping: Sequence[object],
        *,
        multiplicity: object = 1,
        assumptions: object = True,
    ):
        """Validate the parameter map, limits, multiplicity, and assumptions before freezing the region."""
        params = normalize_symbol_sequence(parameters)
        if len(set(params)) != len(params):
            raise ValueError("parametric region parameters must be unique")
        by_name = {param.name: param for param in params}
        if len(by_name) != len(params):
            raise ValueError("parametric region parameter names must be unique")
        sym_limits: list[tuple[sp.Symbol, sp.Expr, sp.Expr]] = []
        seen: set[sp.Symbol] = set()
        for raw, lo, hi in limits:
            if isinstance(raw, str):
                param = by_name.get(raw)
                if param is None:
                    raise ValueError(f"limit variable {raw!r} is not a declared parameter")
            else:
                param = raw
                if param not in params:
                    raise ValueError(f"limit variable {param!r} is not a declared parameter")
            if param in seen:
                raise ValueError(f"duplicate integration limit for parameter {param!r}")
            seen.add(param)
            lower = sp.sympify(lo)
            upper = sp.sympify(hi)
            _validate_interval(lower, upper, label=f"limits for {param}")
            sym_limits.append((param, lower, upper))
        if seen != set(params):
            missing = tuple(param for param in params if param not in seen)
            raise ValueError(f"missing integration limits for parameters {missing!r}")
        mult = sp.sympify(multiplicity)
        try:
            mult_cmp = compare_exact_reals(mult, sp.Integer(0))
        except (TypeError, ValueError, NotImplementedError):
            if mult.is_positive is True:
                mult_cmp = 1
            elif mult.is_nonpositive is True:
                mult_cmp = 0
            else:
                raise ValueError("parametrization multiplicity must be provably positive") from None
        if mult_cmp <= 0:
            raise ValueError("parametrization multiplicity must be positive")
        object.__setattr__(self, "parameters", params)
        object.__setattr__(self, "limits", tuple(sym_limits))
        object.__setattr__(self, "mapping", tuple(sp.sympify(expr) for expr in mapping))
        object.__setattr__(self, "multiplicity", mult)
        object.__setattr__(self, "assumptions", sp.sympify(assumptions))

    def dimension(self) -> int:
        free, _limits, mapping, assumptions = _effective_parametric_data(
            self.parameters, self.limits, self.mapping, self.assumptions
        )
        if assumptions not in (True, sp.true):
            from .regions.operations import region_dimension
            from .symbolic_regions import as_semialgebraic_region

            output = tuple(sp.Dummy(f"image{i + 1}", real=True) for i in range(len(mapping)))
            return region_dimension(as_semialgebraic_region(self, output), output)
        if not free:
            return 0
        return int(sp.Matrix(mapping).jacobian(free).rank())

    def ambient_dimension(self) -> int:
        return len(self.mapping)

    def generic_map_degree(self):
        """Return the generic algebraic fiber degree of this parametrization."""

        from .map_degree import parametric_map_degree

        return parametric_map_degree(self.mapping, self.parameters)


@dataclass(frozen=True)
class TransformedRegion(StandardRegion):
    base: StandardRegion
    mapping: tuple[sp.Expr, ...]
    base_variables: tuple[sp.Symbol, ...]

    def __init__(
        self,
        base: StandardRegion,
        mapping: Sequence[object],
        base_variables: Sequence[sp.Symbol | str],
    ):
        if not isinstance(base, StandardRegion):
            raise TypeError("TransformedRegion base must be a StandardRegion")
        variables = normalize_symbol_sequence(base_variables)
        if len(variables) != base.ambient_dimension():
            raise ValueError("base variable count must match the base ambient dimension")
        object.__setattr__(self, "base", base)
        object.__setattr__(self, "mapping", tuple(sp.sympify(e) for e in mapping))
        object.__setattr__(self, "base_variables", variables)

    def dimension(self) -> int:
        """Compute the exact image dimension from the base tangent space and mapping Jacobian."""
        base_dim = self.base.dimension()
        if base_dim <= 0:
            return max(base_dim, 0)

        substitutions: dict[sp.Symbol, sp.Expr] = {}
        tangent_vectors: tuple[Point, ...] | None = None
        if isinstance(self.base, IntervalRegion):
            if self.base.dimension() == 0:
                return 0
            tangent_vectors = ((sp.Integer(1),),)
        elif isinstance(self.base, BoxRegion):
            vectors: list[Point] = []
            for index, (variable, (lower, upper)) in enumerate(
                zip(self.base_variables, self.base.bounds, strict=True)
            ):
                if _provably_zero(upper - lower):
                    substitutions[variable] = lower
                else:
                    vector = tuple(
                        sp.Integer(1 if i == index else 0) for i in range(len(self.base_variables))
                    )
                    vectors.append(vector)
            tangent_vectors = tuple(vectors)
        elif isinstance(self.base, SimplexRegion):
            anchor = self.base.vertices[0]
            substitutions = dict(zip(self.base_variables, anchor, strict=True))
            differences = tuple(
                tuple(sp.Matrix(vertex) - sp.Matrix(anchor)) for vertex in self.base.vertices[1:]
            )
            tangent_vectors = _independent_vectors(differences)
        elif isinstance(self.base, (ParallelogramRegion, ParallelepipedRegion)):
            substitutions = dict(zip(self.base_variables, self.base.origin, strict=True))
            tangent_vectors = _independent_vectors(self.base.vectors)
        elif isinstance(self.base, PointRegion):
            return 0

        if tangent_vectors is not None:
            if not tangent_vectors:
                return 0
            params = tuple(sp.Dummy(f"u{i + 1}", real=True) for i in range(len(tangent_vectors)))
            base_point = tuple(substitutions.get(v, v) for v in self.base_variables)
            restricted = []
            for expr in self.mapping:
                replacement = {}
                for coord, variable in enumerate(self.base_variables):
                    value = base_point[coord] + sum(
                        param * vector[coord]
                        for param, vector in zip(params, tangent_vectors, strict=True)
                    )
                    replacement[variable] = value
                restricted.append(sp.simplify(expr.subs(replacement)))
            return int(sp.Matrix(restricted).jacobian(params).rank())

        from .regions.operations import region_dimension
        from .symbolic_regions import as_semialgebraic_region

        output = tuple(sp.Dummy(f"image{i + 1}", real=True) for i in range(len(self.mapping)))
        return region_dimension(as_semialgebraic_region(self, output), output)

    def ambient_dimension(self) -> int:
        return len(self.mapping)


@dataclass(frozen=True)
class BooleanRegion(StandardRegion):
    op: str
    regions: tuple[StandardRegion, ...]
    assume_disjoint: bool = False

    def __init__(
        self, op: str, regions: Sequence[StandardRegion], *, assume_disjoint: bool = False
    ):
        if op not in {"union", "intersection", "difference", "symmetric_difference", "complement"}:
            raise ValueError("unsupported BooleanRegion op")
        region_tuple = tuple(regions)
        if any(not isinstance(region, StandardRegion) for region in region_tuple):
            raise TypeError("BooleanRegion members must be StandardRegion objects")
        if op in {"difference", "symmetric_difference"} and len(region_tuple) != 2:
            raise ValueError(f"{op} requires exactly two regions")
        if op == "complement" and len(region_tuple) != 1:
            raise ValueError("complement requires exactly one region")
        ambient_dims = {region.ambient_dimension() for region in region_tuple}
        if len(ambient_dims) > 1:
            raise ValueError("Boolean regions must share one ambient dimension")
        object.__setattr__(self, "op", op)
        object.__setattr__(self, "regions", region_tuple)
        object.__setattr__(self, "assume_disjoint", assume_disjoint)

    def dimension(self) -> int:
        if not self.regions:
            return -1
        if self.op == "union":
            return max(r.dimension() for r in self.regions)
        if self.op == "intersection" and all(
            isinstance(region, IntervalRegion) for region in self.regions
        ):
            intervals = tuple(
                region for region in self.regions if isinstance(region, IntervalRegion)
            )
            lower = intervals[0].lower
            upper = intervals[0].upper
            try:
                for interval in intervals[1:]:
                    if compare_exact_reals(interval.lower, lower) > 0:
                        lower = interval.lower
                    if compare_exact_reals(interval.upper, upper) < 0:
                        upper = interval.upper
                relation = compare_exact_reals(lower, upper)
            except (TypeError, ValueError, NotImplementedError) as exc:
                raise NotImplementedError(
                    "interval-intersection dimension requires exactly comparable endpoints"
                ) from exc
            if relation > 0:
                return -1
            if relation < 0:
                return 1
            for interval in intervals:
                try:
                    at_lower = compare_exact_reals(lower, interval.lower) == 0
                    at_upper = compare_exact_reals(upper, interval.upper) == 0
                except (TypeError, ValueError, NotImplementedError) as exc:
                    raise NotImplementedError(
                        "interval-intersection dimension requires exactly comparable endpoints"
                    ) from exc
                if (at_lower and not interval.lower_closed) or (
                    at_upper and not interval.upper_closed
                ):
                    return -1
            return 0
        raise NotImplementedError(
            "exact dimension is not structurally determined for this BooleanRegion operation; "
            "convert it to SemialgebraicRegion and use region_dimension()"
        )

    def ambient_dimension(self) -> int:
        return self.regions[0].ambient_dimension() if self.regions else 0


def RegionUnion(*regions: StandardRegion, assume_disjoint: bool = False) -> BooleanRegion:
    """Return the Boolean union of standard regions."""
    return BooleanRegion("union", regions, assume_disjoint=assume_disjoint)


def RegionIntersection(*regions: StandardRegion) -> BooleanRegion:
    """Return the Boolean intersection of standard regions."""
    return BooleanRegion("intersection", regions)


def RegionDifference(a: StandardRegion, b: StandardRegion) -> BooleanRegion:
    """Return the Boolean difference of two standard regions."""
    return BooleanRegion("difference", (a, b))


def RegionSymmetricDifference(a: StandardRegion, b: StandardRegion) -> BooleanRegion:
    """Return the Boolean symmetric difference of two standard regions."""
    return BooleanRegion("symmetric_difference", (a, b))


def is_standard_region(obj: object) -> bool:
    return isinstance(obj, StandardRegion)


__all__ = [
    "StandardRegion",
    "PointRegion",
    "IntervalRegion",
    "BoxRegion",
    "SimplexRegion",
    "PolygonRegion",
    "TetrahedronRegion",
    "PolyhedronRegion",
    "ParallelogramRegion",
    "ParallelepipedRegion",
    "PrismRegion",
    "PyramidRegion",
    "BallRegion",
    "SphereRegion",
    "SphericalShellRegion",
    "CylinderRegion",
    "ConeRegion",
    "StadiumRegion",
    "CapsuleRegion",
    "ParametricRegion",
    "TransformedRegion",
    "BooleanRegion",
    "RegionUnion",
    "RegionIntersection",
    "RegionDifference",
    "RegionSymmetricDifference",
    "is_standard_region",
]
