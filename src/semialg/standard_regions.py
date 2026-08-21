from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .exact_arithmetic import compare_exact_reals

Point = tuple[sp.Expr, ...]


def _sympify_point(point: Sequence[object]) -> Point:
    return tuple(sp.sympify(v) for v in point)


def _as_symbols(names: Sequence[sp.Symbol | str]) -> tuple[sp.Symbol, ...]:
    return tuple(sp.Symbol(v, real=True) if isinstance(v, str) else v for v in names)


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


class StandardRegion:
    """Base class for explicit region objects supported by semialg."""

    def dimension(self) -> int:
        raise NotImplementedError

    def ambient_dimension(self) -> int:
        raise NotImplementedError


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
        return 0

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
        return 1 if sp.simplify(self.upper - self.lower) != 0 else 0

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
        return len(self.bounds)

    def ambient_dimension(self) -> int:
        return len(self.bounds)


@dataclass(frozen=True)
class SimplexRegion(StandardRegion):
    vertices: tuple[Point, ...]

    def __init__(self, vertices: Sequence[Sequence[object]]):
        verts = tuple(_sympify_point(v) for v in vertices)
        _validate_same_dimension(verts, label="simplex vertices")
        object.__setattr__(self, "vertices", verts)

    def dimension(self) -> int:
        return max(0, len(self.vertices) - 1)

    def ambient_dimension(self) -> int:
        return len(self.vertices[0]) if self.vertices else 0


@dataclass(frozen=True)
class PolygonRegion(StandardRegion):
    vertices: tuple[Point, ...]

    def __init__(self, vertices: Sequence[Sequence[object]]):
        verts = tuple(_sympify_point(v) for v in vertices)
        if len(verts) < 3:
            raise ValueError("a polygon requires at least three vertices")
        if len({len(v) for v in verts}) != 1 or len(verts[0]) != 2:
            raise ValueError("PolygonRegion currently represents 2D polygons")
        object.__setattr__(self, "vertices", verts)

    def dimension(self) -> int:
        return 2

    def ambient_dimension(self) -> int:
        return 2

    def triangulation(self) -> tuple[SimplexRegion, ...]:
        p0 = self.vertices[0]
        return tuple(
            SimplexRegion((p0, self.vertices[i], self.vertices[i + 1]))
            for i in range(1, len(self.vertices) - 1)
        )


@dataclass(frozen=True)
class TetrahedronRegion(SimplexRegion):
    def __init__(self, vertices: Sequence[Sequence[object]]):
        if len(vertices) != 4:
            raise ValueError("a tetrahedron requires four vertices")
        super().__init__(vertices)


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
        return 3

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
        return 2

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
        return len(self.vectors)

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
        return self.base.dimension() + 1

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
        return self.base.dimension() + 1

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
        return len(self.center)

    def ambient_dimension(self) -> int:
        return len(self.center)


@dataclass(frozen=True)
class SphereRegion(BallRegion):
    def __init__(self, center: Sequence[object], radius: object = 1):
        super().__init__(center, radius)

    def dimension(self) -> int:
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
        return len(self.start)

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
        radius_expr = sp.sympify(radius)
        _validate_nonnegative(radius_expr, label="radius")
        object.__setattr__(self, "start", start_pt)
        object.__setattr__(self, "end", end_pt)
        object.__setattr__(self, "radius", radius_expr)

    def dimension(self) -> int:
        return 2

    def ambient_dimension(self) -> int:
        return 2


@dataclass(frozen=True)
class CapsuleRegion(StadiumRegion):
    def __init__(self, start: Sequence[object], end: Sequence[object], radius: object = 1):
        super().__init__(start, end, radius)

    def dimension(self) -> int:
        return len(self.start)

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
        params = _as_symbols(parameters)
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
        return len(self.parameters)

    def ambient_dimension(self) -> int:
        return len(self.mapping)


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
        object.__setattr__(self, "base", base)
        object.__setattr__(self, "mapping", tuple(sp.sympify(e) for e in mapping))
        object.__setattr__(self, "base_variables", _as_symbols(base_variables))

    def dimension(self) -> int:
        return self.base.dimension()

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
        object.__setattr__(self, "op", op)
        object.__setattr__(self, "regions", tuple(regions))
        object.__setattr__(self, "assume_disjoint", assume_disjoint)

    def dimension(self) -> int:
        if not self.regions:
            return -1
        if self.op == "intersection":
            return min(r.dimension() for r in self.regions)
        return max(r.dimension() for r in self.regions)

    def ambient_dimension(self) -> int:
        return self.regions[0].ambient_dimension() if self.regions else 0


def RegionUnion(*regions: StandardRegion, assume_disjoint: bool = False) -> BooleanRegion:
    return BooleanRegion("union", regions, assume_disjoint=assume_disjoint)


def RegionIntersection(*regions: StandardRegion) -> BooleanRegion:
    return BooleanRegion("intersection", regions)


def RegionDifference(a: StandardRegion, b: StandardRegion) -> BooleanRegion:
    return BooleanRegion("difference", (a, b))


def RegionSymmetricDifference(a: StandardRegion, b: StandardRegion) -> BooleanRegion:
    return BooleanRegion("symmetric_difference", (a, b))


def _volume_unit_ball(n: int) -> sp.Expr:
    return sp.pi ** sp.Rational(n, 2) / sp.gamma(sp.Rational(n, 2) + 1)


def _surface_unit_sphere(n: int) -> sp.Expr:
    return 2 * sp.pi ** sp.Rational(n, 2) / sp.gamma(sp.Rational(n, 2))


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
