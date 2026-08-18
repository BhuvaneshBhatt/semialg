from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

Point = tuple[sp.Expr, ...]


def _sympify_point(point: Sequence[object]) -> Point:
    return tuple(sp.sympify(v) for v in point)


def _as_symbols(names: Sequence[sp.Symbol | str]) -> tuple[sp.Symbol, ...]:
    return tuple(sp.Symbol(v, real=True) if isinstance(v, str) else v for v in names)


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
        object.__setattr__(self, "lower", sp.sympify(lower))
        object.__setattr__(self, "upper", sp.sympify(upper))
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
        object.__setattr__(self, "bounds", tuple((sp.sympify(a), sp.sympify(b)) for a, b in bounds))

    def dimension(self) -> int:
        return len(self.bounds)

    def ambient_dimension(self) -> int:
        return len(self.bounds)


@dataclass(frozen=True)
class SimplexRegion(StandardRegion):
    vertices: tuple[Point, ...]

    def __init__(self, vertices: Sequence[Sequence[object]]):
        object.__setattr__(self, "vertices", tuple(_sympify_point(v) for v in vertices))

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
        object.__setattr__(self, "origin", _sympify_point(origin))
        object.__setattr__(self, "vectors", tuple(_sympify_point(v) for v in vectors))

    def dimension(self) -> int:
        return 2

    def ambient_dimension(self) -> int:
        return len(self.origin)


@dataclass(frozen=True)
class ParallelepipedRegion(StandardRegion):
    origin: Point
    vectors: tuple[Point, ...]

    def __init__(self, origin: Sequence[object], vectors: Sequence[Sequence[object]]):
        object.__setattr__(self, "origin", _sympify_point(origin))
        object.__setattr__(self, "vectors", tuple(_sympify_point(v) for v in vectors))

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
        object.__setattr__(self, "base", base_obj)
        object.__setattr__(self, "vector", _sympify_point(vector))

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
        object.__setattr__(self, "base", base_obj)
        object.__setattr__(self, "apex", _sympify_point(apex))

    def dimension(self) -> int:
        return self.base.dimension() + 1

    def ambient_dimension(self) -> int:
        return len(self.apex)


@dataclass(frozen=True)
class BallRegion(StandardRegion):
    center: Point
    radius: sp.Expr

    def __init__(self, center: Sequence[object], radius: object = 1):
        object.__setattr__(self, "center", _sympify_point(center))
        object.__setattr__(self, "radius", sp.sympify(radius))

    def dimension(self) -> int:
        return len(self.center)

    def ambient_dimension(self) -> int:
        return len(self.center)


@dataclass(frozen=True)
class SphereRegion(BallRegion):
    def dimension(self) -> int:
        return max(0, len(self.center) - 1)


@dataclass(frozen=True)
class SphericalShellRegion(StandardRegion):
    center: Point
    inner_radius: sp.Expr
    outer_radius: sp.Expr

    def __init__(self, center: Sequence[object], radii: tuple[object, object]):
        object.__setattr__(self, "center", _sympify_point(center))
        object.__setattr__(self, "inner_radius", sp.sympify(radii[0]))
        object.__setattr__(self, "outer_radius", sp.sympify(radii[1]))

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
        object.__setattr__(self, "start", _sympify_point(start))
        object.__setattr__(self, "end", _sympify_point(end))
        object.__setattr__(self, "radius", sp.sympify(radius))

    def dimension(self) -> int:
        return len(self.start)

    def ambient_dimension(self) -> int:
        return len(self.start)


@dataclass(frozen=True)
class ConeRegion(CylinderRegion):
    pass


@dataclass(frozen=True)
class StadiumRegion(StandardRegion):
    start: Point
    end: Point
    radius: sp.Expr = sp.Integer(1)

    def __init__(self, start: Sequence[object], end: Sequence[object], radius: object = 1):
        object.__setattr__(self, "start", _sympify_point(start))
        object.__setattr__(self, "end", _sympify_point(end))
        object.__setattr__(self, "radius", sp.sympify(radius))

    def dimension(self) -> int:
        return 2

    def ambient_dimension(self) -> int:
        return 2


@dataclass(frozen=True)
class CapsuleRegion(StadiumRegion):
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
        sym_limits = []
        by_name = {p.name: p for p in params}
        for raw, lo, hi in limits:
            p = by_name.get(raw, sp.Symbol(raw, real=True)) if isinstance(raw, str) else raw
            sym_limits.append((p, sp.sympify(lo), sp.sympify(hi)))
        object.__setattr__(self, "parameters", params)
        object.__setattr__(self, "limits", tuple(sym_limits))
        object.__setattr__(self, "mapping", tuple(sp.sympify(e) for e in mapping))
        object.__setattr__(self, "multiplicity", sp.sympify(multiplicity))
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
