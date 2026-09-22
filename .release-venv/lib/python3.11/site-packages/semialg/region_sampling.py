from __future__ import annotations

import math
import random
from collections.abc import Mapping, Sequence
from functools import lru_cache

import sympy as sp


def _numeric(value: object, *, label: str) -> float:
    expr = sp.sympify(value)
    if expr.free_symbols:
        raise ValueError(f"{label} must be numeric for random sampling")
    result = float(sp.N(expr, 30))
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite for random sampling")
    return result


def _tuple_point(values: Sequence[float], precision: int) -> tuple[sp.Float, ...]:
    return tuple(sp.Float(value, precision) for value in values)


def _unit_direction(dimension: int, rng: random.Random) -> list[float]:
    while True:
        values = [rng.gauss(0.0, 1.0) for _ in range(dimension)]
        norm = math.sqrt(sum(value * value for value in values))
        if norm > 0.0:
            return [value / norm for value in values]


@lru_cache(maxsize=128)
def _ellipsoid_data(region, precision: int):
    center = tuple(_numeric(value, label="ellipsoid center") for value in region.center)
    shape = sp.Matrix(region.shape_matrix).evalf(precision)
    if any(entry.free_symbols for entry in shape):
        raise ValueError("ellipsoid shape matrix must be numeric for random sampling")
    try:
        factor_matrix = shape.cholesky()
    except (ValueError, sp.NonPositiveDefiniteMatrixError) as exc:
        raise ValueError("ellipsoid shape matrix must be positive definite") from exc
    factor = tuple(tuple(float(value) for value in factor_matrix.row(i)) for i in range(shape.rows))
    inverse = factor_matrix.inv().T
    inverse_transpose = tuple(
        tuple(float(value) for value in inverse.row(i)) for i in range(shape.rows)
    )
    envelope = math.sqrt(sum(value * value for row in inverse_transpose for value in row))
    return center, factor, inverse_transpose, envelope


def _matvec(matrix, vector):
    return [
        sum(value * coordinate for value, coordinate in zip(row, vector, strict=True))
        for row in matrix
    ]


def _ellipsoid_boundary_point(region, rng: random.Random, precision: int):
    center, factor, inverse_transpose, envelope = _ellipsoid_data(region, precision)
    # Under x = L u, sphere surface area scales by
    # |det(L)| * ||L^{-T} u||. Rejection sampling with a Frobenius-norm
    # envelope converts a uniform sphere direction into uniform intrinsic
    # surface measure on the ellipsoid.
    dimension = len(center)
    for _ in range(10000):
        direction = _unit_direction(dimension, rng)
        transformed_normal = _matvec(inverse_transpose, direction)
        normal_scale = math.sqrt(sum(value * value for value in transformed_normal))
        if rng.random() * envelope <= normal_scale:
            displacement = _matvec(factor, direction)
            return _tuple_point(
                [origin + delta for origin, delta in zip(center, displacement, strict=True)],
                precision,
            )
    raise ValueError("ellipsoid surface rejection sampling exceeded 10000 attempts")


def _simplex_point(vertices, rng: random.Random, precision: int):
    count = len(vertices)
    if count == 1:
        return tuple(sp.N(value, precision) for value in vertices[0])
    weights = [-math.log(max(rng.random(), 1e-300)) for _ in range(count)]
    total = sum(weights)
    bary = [weight / total for weight in weights]
    ambient = len(vertices[0])
    coords = [
        sum(bary[j] * _numeric(vertices[j][i], label="simplex vertex") for j in range(count))
        for i in range(ambient)
    ]
    return _tuple_point(coords, precision)


@lru_cache(maxsize=128)
def _polytope_sampling_data(region):
    vertices = tuple(region.vertices)
    dimension = region.ambient_dimension()
    lowers = tuple(
        min(_numeric(vertex[i], label="polytope vertex") for vertex in vertices)
        for i in range(dimension)
    )
    uppers = tuple(
        max(_numeric(vertex[i], label="polytope vertex") for vertex in vertices)
        for i in range(dimension)
    )
    hrep = region.h_representation()
    matrix = tuple(
        tuple(float(sp.N(hrep.matrix[row, col], 30)) for col in range(dimension))
        for row in range(hrep.matrix.rows)
    )
    offsets = tuple(float(sp.N(hrep.offsets[row, 0], 30)) for row in range(hrep.matrix.rows))
    return lowers, uppers, matrix, offsets


def _inside_halfspaces(point, matrix, offsets):
    return all(
        sum(coefficient * coordinate for coefficient, coordinate in zip(row, point, strict=True))
        <= offset + 1e-12
        for row, offset in zip(matrix, offsets, strict=True)
    )


@lru_cache(maxsize=128)
def _weighted_region_pieces(region):
    from .standard_regions import Polygon, TetrahedralComplex

    if isinstance(region, Polygon):
        pieces = tuple(region.triangulation())
    elif isinstance(region, TetrahedralComplex):
        pieces = tuple(region.tetrahedra)
    else:
        raise TypeError("weighted decomposition requires a polygon or polyhedron")
    weights = tuple(float(sp.N(piece.measure(), 30)) for piece in pieces)
    return pieces, weights


def _weighted_choice(items, weights, rng: random.Random):
    total = sum(weights)
    if not math.isfinite(total) or total <= 0:
        raise ValueError("sampling weights must have positive finite total")
    target = rng.random() * total
    partial = 0.0
    for item, weight in zip(items, weights, strict=True):
        partial += weight
        if target <= partial:
            return item
    return items[-1]


def _random_standard_region(region, rng: random.Random, precision: int):
    from .standard_regions import (
        Ball,
        Box,
        Ellipsoid,
        EllipsoidBoundary,
        FinitePointSet,
        Interval,
        Parallelepiped,
        Parallelogram,
        Polygon,
        Polytope,
        Simplex,
        Sphere,
        SphericalShell,
        TetrahedralComplex,
    )

    if isinstance(region, FinitePointSet):
        if not region.points:
            raise ValueError("cannot sample an empty point region")
        point = rng.choice(region.points)
        return tuple(sp.N(value, precision) for value in point)

    if isinstance(region, Interval):
        lo = _numeric(region.lower, label="interval lower bound")
        hi = _numeric(region.upper, label="interval upper bound")
        return _tuple_point((lo + (hi - lo) * rng.random(),), precision)

    if isinstance(region, Box):
        values = []
        for index, (lower, upper) in enumerate(region.bounds):
            lo = _numeric(lower, label=f"box lower bound {index}")
            hi = _numeric(upper, label=f"box upper bound {index}")
            values.append(lo + (hi - lo) * rng.random())
        return _tuple_point(values, precision)

    if isinstance(region, Sphere):
        center = [_numeric(value, label="sphere center") for value in region.center]
        radius = _numeric(region.radius, label="sphere radius")
        if radius == 0:
            return _tuple_point(center, precision)
        direction = _unit_direction(len(center), rng)
        return _tuple_point(
            [c + radius * d for c, d in zip(center, direction, strict=True)], precision
        )

    if isinstance(region, Ball):
        center = [_numeric(value, label="ball center") for value in region.center]
        radius = _numeric(region.radius, label="ball radius")
        if radius == 0:
            return _tuple_point(center, precision)
        direction = _unit_direction(len(center), rng)
        radial = radius * rng.random() ** (1.0 / len(center))
        return _tuple_point(
            [c + radial * d for c, d in zip(center, direction, strict=True)], precision
        )

    if isinstance(region, SphericalShell):
        center = [_numeric(value, label="shell center") for value in region.center]
        inner = _numeric(region.inner_radius, label="shell inner radius")
        outer = _numeric(region.outer_radius, label="shell outer radius")
        n = len(center)
        direction = _unit_direction(n, rng)
        radial = (inner**n + rng.random() * (outer**n - inner**n)) ** (1.0 / n)
        return _tuple_point(
            [c + radial * d for c, d in zip(center, direction, strict=True)], precision
        )

    if isinstance(region, EllipsoidBoundary):
        return _ellipsoid_boundary_point(region, rng, precision)

    if isinstance(region, Ellipsoid):
        center, factor, _, _ = _ellipsoid_data(region, precision)
        n = len(center)
        direction = _unit_direction(n, rng)
        radial = rng.random() ** (1.0 / n)
        displacement = _matvec(factor, direction)
        return _tuple_point(
            [origin + radial * delta for origin, delta in zip(center, displacement, strict=True)],
            precision,
        )

    if isinstance(region, Simplex):
        return _simplex_point(region.vertices, rng, precision)

    if isinstance(region, (Polygon, TetrahedralComplex)):
        pieces, weights = _weighted_region_pieces(region)
        piece = _weighted_choice(pieces, weights, rng)
        return _simplex_point(piece.vertices, rng, precision)

    if isinstance(region, Polytope):
        if region.dimension() != region.ambient_dimension():
            raise NotImplementedError(
                "uniform sampling of lower-dimensional Polytope objects is not yet implemented"
            )
        lowers, uppers, matrix, offsets = _polytope_sampling_data(region)
        for _ in range(10000):
            coords = [
                lower + (upper - lower) * rng.random()
                for lower, upper in zip(lowers, uppers, strict=True)
            ]
            if _inside_halfspaces(coords, matrix, offsets):
                return _tuple_point(coords, precision)
        raise ValueError("polytope rejection sampling did not find a point within 10000 attempts")

    if isinstance(region, (Parallelogram, Parallelepiped)):
        if region.dimension() != len(region.vectors):
            raise NotImplementedError(
                "uniform intrinsic sampling requires independent spanning vectors"
            )
        origin = [_numeric(value, label="parallelotope origin") for value in region.origin]
        coords = list(origin)
        for vector in region.vectors:
            coefficient = rng.random()
            for i, value in enumerate(vector):
                coords[i] += coefficient * _numeric(value, label="parallelotope vector")
        return _tuple_point(coords, precision)

    raise NotImplementedError(
        f"uniform random sampling is not implemented for {type(region).__name__}"
    )


def random_points(
    region: object,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    count: int = 1,
    seed: int | None = None,
    distribution: str = "uniform",
    bounds: Sequence[tuple[object, object]]
    | Mapping[sp.Symbol, tuple[object, object]]
    | None = None,
    precision: int = 30,
    attempts: int | None = None,
) -> tuple[tuple[sp.Expr, ...], ...]:
    """Return random points in a region.

    For supported bounded canonical regions, ``distribution="uniform"`` means
    uniform with respect to the region's intrinsic Euclidean (Hausdorff)
    measure. Formula regions use bounded rejection sampling, which is uniform
    with respect to ambient measure on the supplied sampling box.
    """

    if count <= 0:
        return ()
    if precision < 1:
        raise ValueError("precision must be positive")
    key = distribution.lower().replace("-", "_")
    if key != "uniform":
        raise ValueError("only distribution='uniform' is currently supported")

    from .standard_regions import StandardRegion

    rng = random.Random(seed)
    if isinstance(region, StandardRegion):
        return tuple(_random_standard_region(region, rng, precision) for _ in range(count))

    from .sampling import sample_points
    from .symbolic_regions import as_semialgebraic_region

    symbolic = as_semialgebraic_region(region, variables)
    sampled = sample_points(
        symbolic.quantifier_free_formula(),
        symbolic.variables,
        count=count,
        strategy="random",
        exact=False,
        seed=seed,
        bounds=bounds,
        random_attempts=attempts,
        numeric_precision=precision,
    )
    if len(sampled) != count:
        raise ValueError(
            f"random rejection sampling found only {len(sampled)} of {count} requested points"
        )
    return tuple(tuple(point[var] for var in symbolic.variables) for point in sampled)


def random_point(
    region: object,
    variables: Sequence[sp.Symbol | str] | None = None,
    **kwargs,
) -> tuple[sp.Expr, ...]:
    """Return one random point in a region; see :func:`random_points`."""

    points = random_points(region, variables, count=1, **kwargs)
    return points[0]


__all__ = ["random_point", "random_points"]
