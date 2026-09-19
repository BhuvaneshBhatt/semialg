"""Shared exact geometry helpers for standard regions.

Kept separate from the public region classes so validation and elementary
geometry machinery has one owner and the main region module stays navigable.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from ._zero_testing import certified_equal
from .exact_arithmetic import compare_exact_reals

_PointData = tuple[sp.Expr, ...]


def _sympify_point(point: Sequence[object]) -> _PointData:
    return tuple(sp.sympify(v) for v in point)


def _validate_same_dimension(points: Sequence[_PointData], *, label: str) -> int:
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


def _provably_nonpositive(value: sp.Expr) -> bool:
    """Return whether an exact symbolic quantity can be certified as nonpositive."""

    simplified = sp.simplify(value)
    if simplified.is_nonpositive is True:
        return True
    try:
        return compare_exact_reals(simplified, sp.Integer(0)) <= 0
    except (TypeError, ValueError, NotImplementedError):
        return False


def _validate_triangle_sides(a: sp.Expr, b: sp.Expr, c: sp.Expr) -> None:
    """Reject side data that can be certified not to define a nondegenerate triangle."""

    sides = (a, b, c)
    if any(_provably_nonpositive(side) for side in sides):
        raise ValueError("triangle side lengths must be positive")
    for lhs, rhs1, rhs2 in ((a, b, c), (b, a, c), (c, a, b)):
        gap = sp.simplify(rhs1 + rhs2 - lhs)
        if _provably_nonpositive(gap):
            raise ValueError("triangle side lengths must satisfy the strict triangle inequalities")


def _validate_triangle_angle(angle: sp.Expr) -> None:
    """Reject angles that can be certified outside the open interval ``(0, pi)``."""

    if _provably_nonpositive(angle) or _provably_nonpositive(sp.pi - angle):
        raise ValueError("triangle angles must lie strictly between 0 and pi")


def _vector_rank(vectors: Sequence[_PointData]) -> int:
    """Return the exact/generic rank of a finite family of coordinate vectors."""

    if not vectors:
        return 0
    return int(sp.Matrix.hstack(*(sp.Matrix(vector) for vector in vectors)).rank())


def _affine_dimension(points: Sequence[_PointData]) -> int:
    """Return the affine dimension of a finite point set."""

    if len(points) < 2:
        return 0
    anchor = sp.Matrix(points[0])
    vectors = tuple(tuple(sp.Matrix(point) - anchor) for point in points[1:])
    return _vector_rank(vectors)


def _independent_vectors(vectors: Sequence[_PointData]) -> tuple[_PointData, ...]:
    """Return a deterministic maximal independent subfamily."""

    basis: list[_PointData] = []
    rank = 0
    for vector in vectors:
        candidate = (*basis, vector)
        new_rank = _vector_rank(candidate)
        if new_rank > rank:
            basis.append(vector)
            rank = new_rank
    return tuple(basis)


def _identity_directions(dimension: int) -> tuple[_PointData, ...]:
    """Return the standard basis of the ambient coordinate space."""

    return tuple(tuple(sp.Integer(i == j) for i in range(dimension)) for j in range(dimension))


def _sympify_square_matrix(matrix: object, *, dimension: int, label: str) -> sp.ImmutableMatrix:
    """Validate and freeze a square matrix with the requested dimension."""

    value = sp.ImmutableMatrix(matrix)
    if value.shape != (dimension, dimension):
        raise ValueError(f"{label} must be a {dimension} x {dimension} matrix")
    return value


def _positive_definite_conditions(matrix: sp.ImmutableMatrix) -> tuple[sp.Expr, ...]:
    """Return unresolved Sylvester conditions, rejecting a provably invalid matrix."""

    if matrix != matrix.T:
        if any(
            certified_equal(matrix[i, j], matrix[j, i]) is not True
            for i in range(matrix.rows)
            for j in range(i)
        ):
            raise ValueError("ellipsoid shape matrix must be symmetric")
    symmetric = sp.ImmutableMatrix(
        matrix.rows, matrix.cols, lambda i, j: sp.simplify((matrix[i, j] + matrix[j, i]) / 2)
    )
    status = symmetric.is_positive_definite
    if status is False:
        raise ValueError("ellipsoid shape matrix must be positive definite")
    conditions = []
    for size in range(1, symmetric.rows + 1):
        condition = sp.simplify(symmetric[:size, :size].det() > 0)
        if condition is sp.false:
            raise ValueError("ellipsoid shape matrix must be positive definite")
        if condition is not sp.true:
            conditions.append(condition)
    return tuple(conditions)


def _sphere_through_points(points: Sequence[Sequence[object]]) -> tuple[_PointData, sp.Expr]:
    """Return the unique sphere center/radius through full-dimensional points."""

    pts = tuple(_sympify_point(point) for point in points)
    if not pts:
        raise ValueError("sphere construction requires points")
    ambient = _validate_same_dimension(pts, label="sphere points")
    if ambient == 0:
        raise ValueError("sphere points must have positive ambient dimension")
    if len(pts) != ambient + 1:
        raise ValueError("a unique sphere in n dimensions requires exactly n + 1 points")
    if _affine_dimension(pts) != ambient:
        raise ValueError("sphere points must be affinely independent")
    p0 = sp.Matrix(pts[0])
    rows = []
    rhs = []
    for point in pts[1:]:
        pi = sp.Matrix(point)
        rows.append(tuple(2 * (pi - p0)))
        rhs.append(sp.expand(pi.dot(pi) - p0.dot(p0)))
    matrix = sp.Matrix(rows)
    center_vec = matrix.LUsolve(sp.Matrix(rhs))
    center = tuple(sp.simplify(value) for value in center_vec)
    radius = sp.simplify(sp.sqrt(sum((center[i] - pts[0][i]) ** 2 for i in range(ambient))))
    return center, radius


def _simplex_face_measure(vertices: Sequence[_PointData]) -> sp.Expr:
    """Return the intrinsic Euclidean measure of a simplex face."""

    dim = len(vertices) - 1
    if dim == 0:
        return sp.Integer(1)
    anchor = sp.Matrix(vertices[0])
    columns = [sp.Matrix(vertex) - anchor for vertex in vertices[1:]]
    gram = sp.Matrix.hstack(*columns).T * sp.Matrix.hstack(*columns)
    return sp.simplify(sp.sqrt(gram.det()) / sp.factorial(dim))


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


def _orientation(a: _PointData, b: _PointData, c: _PointData) -> int:
    """Return the certified orientation sign of three exact planar points."""

    cross = sp.expand((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
    try:
        return compare_exact_reals(cross, sp.Integer(0))
    except (TypeError, ValueError, NotImplementedError) as exc:
        raise ValueError("polygon orientation could not be certified exactly") from exc


def _on_segment(a: _PointData, b: _PointData, p: _PointData) -> bool:
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


def _segments_intersect(a: _PointData, b: _PointData, c: _PointData, d: _PointData) -> bool:
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


def _validate_simple_polygon(vertices: Sequence[_PointData]) -> None:
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


def _signed_polygon_orientation(vertices: Sequence[_PointData]) -> int:
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


def _point_in_triangle(
    point: _PointData, triangle: tuple[_PointData, _PointData, _PointData], sign: int
) -> bool:
    a, b, c = triangle
    return all(
        sign * orient >= 0
        for orient in (
            _orientation(a, b, point),
            _orientation(b, c, point),
            _orientation(c, a, point),
        )
    )


def _triangulate_polygon(
    vertices: Sequence[_PointData],
) -> tuple[tuple[_PointData, _PointData, _PointData], ...]:
    """Triangulate a certified simple polygon by exact ear clipping."""

    sign = _signed_polygon_orientation(vertices)
    remaining = list(range(len(vertices)))
    triangles: list[tuple[_PointData, _PointData, _PointData]] = []
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
