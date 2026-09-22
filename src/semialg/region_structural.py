"""Fast exact rules for relations and operations on canonical regions."""

from __future__ import annotations

import sympy as sp

from .exact_arithmetic import compare_exact_reals
from .standard_regions import (
    Ball,
    Box,
    Ellipsoid,
    FinitePointSet,
    Hyperplane,
    Interval,
    Line,
    Parallelepiped,
    Point,
    Polygon,
    Polytope,
    Ray,
    Simplex,
    Sphere,
    StandardRegion,
    TransformedRegion,
)

UNRESOLVED = object()


def _compare(left: sp.Expr, right: sp.Expr) -> int | None:
    try:
        return compare_exact_reals(sp.simplify(left), sp.simplify(right))
    except (TypeError, ValueError, NotImplementedError):
        return None


def _filled_ball(region: object) -> bool:
    return isinstance(region, Ball) and not isinstance(region, Sphere)


def _interval_empty(region: Interval) -> bool:
    relation = _compare(region.lower, region.upper)
    return relation == 0 and not (region.lower_closed and region.upper_closed)


def _point_in_interval(point: sp.Expr, region: Interval) -> bool | None:
    lo = _compare(point, region.lower)
    hi = _compare(point, region.upper)
    if lo is None or hi is None:
        return None
    if lo < 0 or hi > 0:
        return False
    if lo == 0 and not region.lower_closed:
        return False
    if hi == 0 and not region.upper_closed:
        return False
    return True


def _interval_subset(left: Interval, right: Interval) -> bool | None:
    if _interval_empty(left):
        return True
    if _interval_empty(right):
        return _interval_empty(left)
    lo = _compare(left.lower, right.lower)
    hi = _compare(left.upper, right.upper)
    if lo is None or hi is None:
        return None
    if lo < 0 or hi > 0:
        return False
    if lo == 0 and left.lower_closed and not right.lower_closed:
        return False
    if hi == 0 and left.upper_closed and not right.upper_closed:
        return False
    return True


def _interval_disjoint(left: Interval, right: Interval) -> bool | None:
    if _interval_empty(left) or _interval_empty(right):
        return True
    left_before = _compare(left.upper, right.lower)
    if left_before is not None:
        if left_before < 0:
            return True
        if left_before == 0:
            return not (left.upper_closed and right.lower_closed)
    right_before = _compare(right.upper, left.lower)
    if right_before is not None:
        if right_before < 0:
            return True
        if right_before == 0:
            return not (right.upper_closed and left.lower_closed)
    if left_before is not None and right_before is not None:
        return False
    return None


def _polytope_constraints(region: Polytope):
    if region.dimension() != region.ambient_dimension():
        return None
    try:
        return region.h_representation()
    except (TypeError, ValueError, NotImplementedError):
        return None


def _box_hrep(region: Box):
    from .polyhedral import HRepresentation

    n = region.ambient_dimension()
    matrix = []
    offsets = []
    for i, (lower, upper) in enumerate(region.bounds):
        upper_row = [sp.Integer(0)] * n
        upper_row[i] = sp.Integer(1)
        matrix.append(upper_row)
        offsets.append(upper)
        lower_row = [sp.Integer(0)] * n
        lower_row[i] = sp.Integer(-1)
        matrix.append(lower_row)
        offsets.append(-lower)
    return HRepresentation(matrix, offsets)


def _as_polytope(region):
    if isinstance(region, Polytope):
        return region
    if isinstance(region, (Simplex, Polygon, Parallelepiped)):
        try:
            return Polytope.from_region(region)
        except (TypeError, ValueError, NotImplementedError):
            return None
    if isinstance(region, Box):
        import itertools

        vertices = tuple(itertools.product(*((lower, upper) for lower, upper in region.bounds)))
        return Polytope(vertices)
    return None


def _convex_hrep(region):
    if isinstance(region, Box):
        return _box_hrep(region)
    polytope = _as_polytope(region)
    return None if polytope is None else _polytope_constraints(polytope)


def _point_satisfies_hrep(point, representation) -> bool | None:
    vector = sp.Matrix(point)
    for row in range(representation.matrix.rows):
        residual = sp.simplify(
            (representation.matrix.row(row) * vector)[0] - representation.offsets[row, 0]
        )
        relation = _compare(residual, sp.Integer(0))
        if relation is None:
            return None
        if relation > 0:
            return False
    return True


def _clip_parameter_interval(point, direction, representation, *, ray: bool = False):
    lower = sp.Integer(0) if ray else -sp.oo
    upper = sp.oo
    p = sp.Matrix(point)
    d = sp.Matrix(direction)
    for row in range(representation.matrix.rows):
        normal = representation.matrix.row(row)
        alpha = sp.simplify((normal * d)[0])
        beta = sp.simplify(representation.offsets[row, 0] - (normal * p)[0])
        sign = _compare(alpha, sp.Integer(0))
        if sign is None:
            return None
        if sign == 0:
            feasible = _compare(beta, sp.Integer(0))
            if feasible is None:
                return None
            if feasible < 0:
                return False
            continue
        bound = sp.simplify(beta / alpha)
        if sign > 0:
            if upper is sp.oo:
                upper = bound
            else:
                cmp = _compare(bound, upper)
                if cmp is None:
                    return None
                if cmp < 0:
                    upper = bound
        else:
            if lower is -sp.oo:
                lower = bound
            else:
                cmp = _compare(bound, lower)
                if cmp is None:
                    return None
                if cmp > 0:
                    lower = bound
    if lower is not -sp.oo and upper is not sp.oo:
        order = _compare(lower, upper)
        if order is None:
            return None
        if order > 0:
            return False
    return lower, upper


def _parameterized_segment(point, direction, bounds):
    lower, upper = bounds
    t = sp.Dummy("t", real=True)
    mapping = tuple(sp.simplify(a + t * b) for a, b in zip(point, direction, strict=True))
    if lower is -sp.oo or upper is sp.oo:
        return UNRESOLVED
    return TransformedRegion(Interval(lower, upper), mapping, (t,))


def _line_convex_intersection(line, region):
    representation = _convex_hrep(region)
    if representation is None:
        return UNRESOLVED
    bounds = _clip_parameter_interval(
        line.point, line.direction, representation, ray=isinstance(line, Ray)
    )
    if bounds is None or bounds is False:
        return UNRESOLVED
    return _parameterized_segment(line.point, line.direction, bounds)


def _hyperplane_polytope_intersection(plane: Hyperplane, region):
    polytope = _as_polytope(region)
    if polytope is None or polytope.dimension() != polytope.ambient_dimension():
        return UNRESOLVED
    if plane.ambient_dimension() != polytope.ambient_dimension():
        raise ValueError("region ambient dimensions do not match")
    normal = sp.Matrix(plane.normal)
    anchor = sp.Matrix(plane.point)
    values = [sp.simplify(normal.dot(sp.Matrix(v) - anchor)) for v in polytope.vertices]
    signs = [_compare(value, sp.Integer(0)) for value in values]
    if any(sign is None for sign in signs):
        return UNRESOLVED
    if all(sign > 0 for sign in signs) or all(sign < 0 for sign in signs):
        return UNRESOLVED
    points = [v for v, sign in zip(polytope.vertices, signs, strict=True) if sign == 0]
    try:
        adjacency = polytope.adjacency().vertex_neighbors
    except (TypeError, ValueError, NotImplementedError):
        return UNRESOLVED
    seen = set()
    for i, neighbors in enumerate(adjacency):
        for j in neighbors:
            edge = tuple(sorted((i, j)))
            if edge in seen:
                continue
            seen.add(edge)
            si, sj = signs[i], signs[j]
            if si == 0 or sj == 0 or si * sj >= 0:
                continue
            vi = sp.Matrix(polytope.vertices[i])
            vj = sp.Matrix(polytope.vertices[j])
            denominator = sp.simplify(values[i] - values[j])
            t = sp.simplify(values[i] / denominator)
            point = tuple(sp.simplify(x) for x in (vi + t * (vj - vi)))
            if point not in points:
                points.append(point)
    if not points:
        return UNRESOLVED
    return Polytope(points)


def _hyperplane_ball_intersection(plane: Hyperplane, ball):
    if not _filled_ball(ball):
        return UNRESOLVED
    if plane.ambient_dimension() != ball.ambient_dimension():
        raise ValueError("region ambient dimensions do not match")
    normal = sp.Matrix(plane.normal)
    center = sp.Matrix(ball.center)
    anchor = sp.Matrix(plane.point)
    norm_sq = sp.simplify(normal.dot(normal))
    delta = sp.simplify(normal.dot(center - anchor))
    radius_sq = sp.simplify(ball.radius**2 - delta**2 / norm_sq)
    sign = _compare(radius_sq, sp.Integer(0))
    if sign is None or sign < 0:
        return UNRESOLVED
    closest = center - normal * sp.simplify(delta / norm_sq)
    if sign == 0:
        return Point(tuple(sp.simplify(value) for value in closest))
    basis = sp.Matrix([plane.normal]).nullspace()
    if not basis:
        return Point(tuple(sp.simplify(value) for value in closest))
    orthonormal = sp.GramSchmidt(basis, True)
    variables = tuple(sp.Dummy(f"u{i}", real=True) for i in range(len(orthonormal)))
    mapping = tuple(
        sp.simplify(
            closest[row] + sum(v * b[row] for v, b in zip(variables, orthonormal, strict=True))
        )
        for row in range(plane.ambient_dimension())
    )
    base = Ball((0,) * len(variables), sp.sqrt(radius_sq))
    return TransformedRegion(base, mapping, variables)


def _hyperplane_ellipsoid_intersection(plane: Hyperplane, ellipsoid: Ellipsoid):
    if plane.ambient_dimension() != ellipsoid.ambient_dimension():
        raise ValueError("region ambient dimensions do not match")
    normal = sp.Matrix(plane.normal)
    center = sp.Matrix(ellipsoid.center)
    anchor = sp.Matrix(plane.point)
    shape = sp.Matrix(ellipsoid.shape_matrix)
    denominator = sp.simplify((normal.T * shape * normal)[0])
    denominator_sign = _compare(denominator, sp.Integer(0))
    if denominator_sign is None or denominator_sign <= 0:
        return UNRESOLVED
    delta = sp.simplify(normal.dot(center - anchor))
    residual = sp.simplify(1 - delta**2 / denominator)
    sign = _compare(residual, sp.Integer(0))
    if sign is None or sign < 0:
        return UNRESOLVED
    section_center = center - shape * normal * sp.simplify(delta / denominator)
    if sign == 0:
        return Point(tuple(sp.simplify(value) for value in section_center))
    basis = sp.Matrix([plane.normal]).nullspace()
    if not basis:
        return Point(tuple(sp.simplify(value) for value in section_center))
    B = sp.Matrix.hstack(*basis)
    restricted = sp.simplify(B.T * shape.inv() * B)
    section_shape = sp.ImmutableMatrix(sp.simplify(residual) * restricted.inv())
    variables = tuple(sp.Dummy(f"u{i}", real=True) for i in range(B.cols))
    mapping = tuple(
        sp.simplify(
            section_center[row] + sum(B[row, col] * variables[col] for col in range(B.cols))
        )
        for row in range(plane.ambient_dimension())
    )
    base = Ellipsoid((0,) * B.cols, section_shape)
    return TransformedRegion(base, mapping, variables)


def _line_ball_intersection(line, ball):
    if not _filled_ball(ball):
        return UNRESOLVED
    p = sp.Matrix(line.point) - sp.Matrix(ball.center)
    d = sp.Matrix(line.direction)
    a = sp.simplify(d.dot(d))
    b = sp.simplify(2 * p.dot(d))
    c = sp.simplify(p.dot(p) - ball.radius**2)
    discriminant = sp.factor(b**2 - 4 * a * c)
    disc_sign = _compare(discriminant, sp.Integer(0))
    if disc_sign is None or disc_sign < 0:
        return UNRESOLVED
    root = sp.sqrt(discriminant)
    lo = sp.simplify((-b - root) / (2 * a))
    hi = sp.simplify((-b + root) / (2 * a))
    if isinstance(line, Ray):
        hi_sign = _compare(hi, sp.Integer(0))
        if hi_sign is None or hi_sign < 0:
            return UNRESOLVED
        lo_sign = _compare(lo, sp.Integer(0))
        if lo_sign is None:
            return UNRESOLVED
        if lo_sign < 0:
            lo = sp.Integer(0)
    return _parameterized_segment(line.point, line.direction, (lo, hi))


def structural_subset(left: object, right: object) -> bool | None:
    """Return a certified subset result for canonical cases, or ``None``."""

    if left == right:
        return True
    if not isinstance(left, StandardRegion) or not isinstance(right, StandardRegion):
        return None
    if left.ambient_dimension() != right.ambient_dimension():
        raise ValueError("region ambient dimensions do not match")
    if isinstance(left, FinitePointSet):
        try:
            values = tuple(right.contains(point) for point in left.points)
        except (TypeError, ValueError, NotImplementedError):
            return None
        return all(values)
    if isinstance(left, Interval) and isinstance(right, Interval):
        return _interval_subset(left, right)
    if isinstance(left, Box) and isinstance(right, Box):
        if len(left.bounds) != len(right.bounds):
            return False
        results = [
            _interval_subset(Interval(*a), Interval(*b))
            for a, b in zip(left.bounds, right.bounds, strict=True)
        ]
        if any(result is False for result in results):
            return False
        return True if all(result is True for result in results) else None
    if isinstance(left, Polytope) and isinstance(right, Polytope):
        representation = _polytope_constraints(right)
        if representation is not None:
            results = [_point_satisfies_hrep(vertex, representation) for vertex in left.vertices]
            if any(result is False for result in results):
                return False
            if all(result is True for result in results):
                return True
            return None
    if _filled_ball(left) and _filled_ball(right):
        radius_gap = _compare(right.radius, left.radius)
        if radius_gap is not None and radius_gap < 0:
            return False
        if radius_gap is None:
            return None
        distance_sq = sp.Add(
            *((a - b) ** 2 for a, b in zip(left.center, right.center, strict=True))
        )
        allowance = sp.expand((right.radius - left.radius) ** 2)
        relation = _compare(distance_sq, allowance)
        return None if relation is None else relation <= 0
    return None


def structural_disjoint(left: object, right: object) -> bool | None:
    """Return a certified disjointness result for canonical cases, or ``None``."""

    if not isinstance(left, StandardRegion) or not isinstance(right, StandardRegion):
        return None
    if left == right:
        return left.dimension() < 0
    if left.ambient_dimension() != right.ambient_dimension():
        raise ValueError("region ambient dimensions do not match")
    if isinstance(left, FinitePointSet):
        try:
            values = tuple(not right.contains(point) for point in left.points)
        except (TypeError, ValueError, NotImplementedError):
            return None
        return all(values)
    if isinstance(right, FinitePointSet):
        return structural_disjoint(right, left)
    if isinstance(left, Interval) and isinstance(right, Interval):
        return _interval_disjoint(left, right)
    if isinstance(left, Box) and isinstance(right, Box):
        results = [
            _interval_disjoint(Interval(*a), Interval(*b))
            for a, b in zip(left.bounds, right.bounds, strict=True)
        ]
        if any(result is True for result in results):
            return True
        return False if all(result is False for result in results) else None
    if isinstance(left, Polytope) and isinstance(right, Polytope):
        for source, target in ((left, right), (right, left)):
            representation = _polytope_constraints(source)
            if representation is None:
                continue
            for row in range(representation.matrix.rows):
                residuals = [
                    sp.simplify(
                        (representation.matrix.row(row) * sp.Matrix(vertex))[0]
                        - representation.offsets[row, 0]
                    )
                    for vertex in target.vertices
                ]
                signs = [_compare(value, sp.Integer(0)) for value in residuals]
                if all(sign is not None and sign > 0 for sign in signs):
                    return True
        return None
    if _filled_ball(left) and _filled_ball(right):
        distance_sq = sp.Add(
            *((a - b) ** 2 for a, b in zip(left.center, right.center, strict=True))
        )
        separation_sq = sp.expand((left.radius + right.radius) ** 2)
        relation = _compare(distance_sq, separation_sq)
        return None if relation is None else relation > 0
    return None


def structural_intersection(left: object, right: object):
    """Return a canonical intersection when one is cheaply determined."""

    if left == right:
        return left
    convex_polyhedral = (Box, Simplex, Polygon, Parallelepiped, Polytope)
    if isinstance(left, (Line, Ray)) and isinstance(right, convex_polyhedral):
        return _line_convex_intersection(left, right)
    if isinstance(right, (Line, Ray)) and isinstance(left, convex_polyhedral):
        return _line_convex_intersection(right, left)
    if isinstance(left, (Line, Ray)) and _filled_ball(right):
        return _line_ball_intersection(left, right)
    if isinstance(right, (Line, Ray)) and _filled_ball(left):
        return _line_ball_intersection(right, left)
    if isinstance(left, Hyperplane) and isinstance(right, convex_polyhedral):
        return _hyperplane_polytope_intersection(left, right)
    if isinstance(right, Hyperplane) and isinstance(left, convex_polyhedral):
        return _hyperplane_polytope_intersection(right, left)
    if isinstance(left, Hyperplane) and _filled_ball(right):
        return _hyperplane_ball_intersection(left, right)
    if isinstance(right, Hyperplane) and _filled_ball(left):
        return _hyperplane_ball_intersection(right, left)
    if isinstance(left, Hyperplane) and isinstance(right, Ellipsoid):
        return _hyperplane_ellipsoid_intersection(left, right)
    if isinstance(right, Hyperplane) and isinstance(left, Ellipsoid):
        return _hyperplane_ellipsoid_intersection(right, left)
    if isinstance(left, Interval) and isinstance(right, Interval):
        if _interval_empty(left):
            return left
        if _interval_empty(right):
            return right
        if _interval_disjoint(left, right) is True:
            return UNRESOLVED
        lo_cmp = _compare(left.lower, right.lower)
        hi_cmp = _compare(left.upper, right.upper)
        if lo_cmp is None or hi_cmp is None:
            return UNRESOLVED
        lower = left.lower if lo_cmp >= 0 else right.lower
        upper = left.upper if hi_cmp <= 0 else right.upper
        lower_closed = (
            (left.lower_closed if lo_cmp >= 0 else right.lower_closed)
            if lo_cmp != 0
            else left.lower_closed and right.lower_closed
        )
        upper_closed = (
            (left.upper_closed if hi_cmp <= 0 else right.upper_closed)
            if hi_cmp != 0
            else left.upper_closed and right.upper_closed
        )
        if _compare(lower, upper) == 0 and not (lower_closed and upper_closed):
            return UNRESOLVED
        return Interval(lower, upper, lower_closed=lower_closed, upper_closed=upper_closed)
    if isinstance(left, Polytope) and isinstance(right, Polytope):
        try:
            from .polyhedral_booleans import polyhedral_intersection

            result = polyhedral_intersection(left, right)
        except (ValueError, NotImplementedError):
            result = None
        if result is not None:
            return result
    if isinstance(left, Box) and isinstance(right, Box):
        if len(left.bounds) != len(right.bounds):
            return UNRESOLVED
        bounds = []
        for a, b in zip(left.bounds, right.bounds, strict=True):
            piece = structural_intersection(Interval(*a), Interval(*b))
            if not isinstance(piece, Interval):
                return UNRESOLVED
            bounds.append((piece.lower, piece.upper))
        return Box(bounds)
    subset = structural_subset(left, right)
    if subset is True:
        return left
    reverse = structural_subset(right, left)
    if reverse is True:
        return right
    return UNRESOLVED


def structural_product(left: object, right: object):
    """Return a canonical Cartesian product when one has a direct representation."""

    if isinstance(left, Interval) and isinstance(right, Interval):
        if not all((left.lower_closed, left.upper_closed, right.lower_closed, right.upper_closed)):
            return UNRESOLVED
        return Box(((left.lower, left.upper), (right.lower, right.upper)))
    if isinstance(left, Box) and isinstance(right, Box):
        return Box((*left.bounds, *right.bounds))
    if isinstance(left, Interval) and isinstance(right, Box):
        if not (left.lower_closed and left.upper_closed):
            return UNRESOLVED
        return Box(((left.lower, left.upper), *right.bounds))
    if isinstance(left, Box) and isinstance(right, Interval):
        if not (right.lower_closed and right.upper_closed):
            return UNRESOLVED
        return Box((*left.bounds, (right.lower, right.upper)))
    if isinstance(left, Point) and isinstance(right, Point):
        return Point((*left.coordinates, *right.coordinates))
    return UNRESOLVED


def structural_minkowski_sum(left: object, right: object):
    """Return a canonical Minkowski sum when a closed form is available."""

    if isinstance(left, Point) and isinstance(right, Point):
        return Point(tuple(a + b for a, b in zip(left.coordinates, right.coordinates, strict=True)))
    if isinstance(left, Point) and isinstance(right, StandardRegion):
        if left.ambient_dimension() == right.ambient_dimension():
            return right.transform(sp.eye(right.ambient_dimension()), left.coordinates)
    if isinstance(right, Point) and isinstance(left, StandardRegion):
        return structural_minkowski_sum(right, left)
    if isinstance(left, Interval) and isinstance(right, Interval):
        return Interval(
            left.lower + right.lower,
            left.upper + right.upper,
            lower_closed=left.lower_closed and right.lower_closed,
            upper_closed=left.upper_closed and right.upper_closed,
        )
    if isinstance(left, Box) and isinstance(right, Box):
        if len(left.bounds) != len(right.bounds):
            raise ValueError("region ambient dimensions do not match")
        return Box(
            tuple(
                (a0 + b0, a1 + b1)
                for (a0, a1), (b0, b1) in zip(left.bounds, right.bounds, strict=True)
            )
        )
    if isinstance(left, Polytope) and isinstance(right, Polytope):
        if left.ambient_dimension() != right.ambient_dimension():
            raise ValueError("region ambient dimensions do not match")
        sums = tuple(
            tuple(sp.simplify(a + b) for a, b in zip(v, w, strict=True))
            for v in left.vertices
            for w in right.vertices
        )
        return Polytope(sums)
    if _filled_ball(left) and _filled_ball(right):
        return Ball(
            tuple(a + b for a, b in zip(left.center, right.center, strict=True)),
            left.radius + right.radius,
        )
    return UNRESOLVED
