"""Structure-preserving affine and symbolic transformations of canonical regions."""

from __future__ import annotations

import itertools

import sympy as sp

from ._zero_testing import certified_equal, certified_nonzero, certified_zero
from .internal_symbols import collision_free_real_symbols
from .polyhedral import HRepresentation
from .standard_regions import (
    AffineHalfSpace,
    AffineSpace,
    Ball,
    ConicRegion,
    Ellipsoid,
    EllipsoidBoundary,
    HalfSpace,
    Hyperplane,
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


def _require_real_affine_data(*matrices: sp.MatrixBase) -> None:
    """Reject affine coefficients that are not certified real."""
    for matrix in matrices:
        for entry in matrix:
            if sp.ask(sp.Q.real(entry)) is not True:
                raise ValueError("affine map matrix and offset entries must be real-valued")


def _map_data(matrix, offset, input_dimension: int):
    A = sp.ImmutableMatrix(matrix)
    if A.cols != input_dimension:
        raise ValueError("affine map column count must match the region ambient dimension")
    b = sp.ImmutableMatrix.zeros(A.rows, 1) if offset is None else sp.ImmutableMatrix(offset)
    if b.shape != (A.rows, 1):
        raise ValueError("affine offset must match the map output dimension")
    _require_real_affine_data(A, b)
    return A, b


def _mapped_point(A: sp.MatrixBase, b: sp.MatrixBase, point) -> tuple[sp.Expr, ...]:
    value = A * sp.Matrix(point) + b
    return tuple(sp.simplify(value[i, 0]) for i in range(A.rows))


def _mapped_vector(A: sp.MatrixBase, vector) -> tuple[sp.Expr, ...]:
    value = A * sp.Matrix(vector)
    return tuple(sp.simplify(value[i, 0]) for i in range(A.rows))


def _is_zero_vector(vector) -> bool:
    return all(certified_zero(value) is True for value in vector)


def _is_invertible_square(A: sp.MatrixBase) -> bool:
    if A.rows != A.cols:
        return False
    return certified_nonzero(A.det())


def _similarity_scale_squared(A: sp.MatrixBase):
    """Return ``s**2`` when a square map is certifiably a Euclidean similarity."""
    if A.rows != A.cols:
        return None
    gram = sp.ImmutableMatrix(A.T * A)
    scale_sq = sp.simplify(gram[0, 0]) if gram.rows else sp.Integer(1)
    for i in range(gram.rows):
        for j in range(gram.cols):
            expected = scale_sq if i == j else sp.Integer(0)
            if certified_equal(gram[i, j], expected) is not True:
                return None
    if not certified_nonzero(scale_sq):
        return None
    return scale_sq


def _lazy_affine_image(region: StandardRegion, A: sp.MatrixBase, b: sp.MatrixBase):
    variables = sp.symbols(f"u0:{region.ambient_dimension()}", real=True)
    mapping = tuple(A * sp.Matrix(variables) + b)
    return TransformedRegion(region, mapping, variables)


def affine_image(region, matrix, offset=None):
    """Return the exact affine image ``A*x+b``.

    The result stays canonical whenever the image belongs to a canonical family.
    Rank-changing maps and maps whose rank cannot be proved remain exact as
    :class:`TransformedRegion` values rather than assuming invertibility.
    """

    if isinstance(region, HRepresentation):
        A, b = _map_data(matrix, offset, region.ambient_dimension)
        if not _is_invertible_square(A):
            return affine_image(region.to_polytope(), A, tuple(b))
        inv = A.inv()
        new_matrix = region.matrix * inv
        new_offsets = region.offsets + new_matrix * b
        return HRepresentation(new_matrix.tolist(), tuple(new_offsets))

    if not isinstance(region, StandardRegion):
        raise TypeError("affine_image requires a canonical region or HRepresentation")
    A, b = _map_data(matrix, offset, region.ambient_dimension())

    if isinstance(region, Point):
        return Point(_mapped_point(A, b, region.coordinates))
    if isinstance(region, Line):
        point = _mapped_point(A, b, region.point)
        direction = _mapped_vector(A, region.direction)
        return Point(point) if _is_zero_vector(direction) else Line(point, direction)
    if isinstance(region, Ray):
        point = _mapped_point(A, b, region.point)
        direction = _mapped_vector(A, region.direction)
        return Point(point) if _is_zero_vector(direction) else Ray(point, direction)
    if isinstance(region, AffineSpace):
        return AffineSpace(
            _mapped_point(A, b, region.point),
            tuple(_mapped_vector(A, direction) for direction in region.directions),
        )
    if isinstance(region, AffineHalfSpace):
        point = _mapped_point(A, b, region.point)
        directions = tuple(_mapped_vector(A, direction) for direction in region.directions)
        inward = _mapped_vector(A, region.inward)
        if _is_zero_vector(inward):
            return AffineSpace(point, directions)
        try:
            return AffineHalfSpace(point, directions, inward)
        except ValueError:
            return _lazy_affine_image(region, A, b)
    if isinstance(region, (Hyperplane, HalfSpace)):
        if not _is_invertible_square(A):
            return _lazy_affine_image(region, A, b)
        point = _mapped_point(A, b, region.point)
        normal = tuple(A.inv().T * sp.Matrix(region.normal))
        cls = Hyperplane if isinstance(region, Hyperplane) else HalfSpace
        return cls(normal, point)
    if isinstance(region, Simplex):
        vertices = tuple(_mapped_point(A, b, vertex) for vertex in region.vertices)
        try:
            return Simplex(vertices, construction_conditions=region.construction_conditions)
        except ValueError:
            return Polytope(vertices)
    if isinstance(region, Polygon):
        vertices = tuple(_mapped_point(A, b, vertex) for vertex in region.vertices)
        if A.rows == 2 and _is_invertible_square(A):
            return Polygon(vertices)
        return Polytope(vertices)
    if isinstance(region, Parallelepiped):
        origin = _mapped_point(A, b, region.origin)
        vectors = tuple(_mapped_vector(A, vector) for vector in region.vectors)
        try:
            return Parallelepiped(origin, vectors)
        except ValueError:
            vertices = []
            for mask in itertools.product((0, 1), repeat=len(region.vectors)):
                vertices.append(
                    tuple(
                        sp.simplify(
                            origin[j]
                            + sum(bit * v[j] for bit, v in zip(mask, vectors, strict=True))
                        )
                        for j in range(A.rows)
                    )
                )
            return Polytope(vertices)
    if isinstance(region, Polytope):
        return Polytope(tuple(_mapped_point(A, b, vertex) for vertex in region.vertices))
    if isinstance(region, ConicRegion):
        point = _mapped_point(A, b, region.point)
        directions = tuple(_mapped_vector(A, d) for d in region.directions)
        rays = tuple(_mapped_vector(A, r) for r in region.rays)
        directions = tuple(d for d in directions if not _is_zero_vector(d))
        rays = tuple(r for r in rays if not _is_zero_vector(r))
        return ConicRegion(point, directions=directions, rays=rays)

    if isinstance(region, (Ball, Sphere)):
        center = _mapped_point(A, b, region.center)
        if certified_zero(region.radius) is True:
            return Point(center)
        scale_sq = _similarity_scale_squared(A)
        if scale_sq is not None:
            cls = Sphere if isinstance(region, Sphere) else Ball
            return cls(center, sp.simplify(region.radius * sp.sqrt(scale_sq)))
        if _is_invertible_square(A):
            shape = sp.ImmutableMatrix(A * region.shape_matrix * A.T)
            cls = EllipsoidBoundary if isinstance(region, Sphere) else Ellipsoid
            return cls(center, shape)
        return _lazy_affine_image(region, A, b)

    if isinstance(region, (Ellipsoid, EllipsoidBoundary)):
        if _is_invertible_square(A):
            center = _mapped_point(A, b, region.center)
            shape = sp.ImmutableMatrix(A * region.shape_matrix * A.T)
            cls = EllipsoidBoundary if isinstance(region, EllipsoidBoundary) else Ellipsoid
            return cls(center, shape)
        return _lazy_affine_image(region, A, b)

    return _lazy_affine_image(region, A, b)


def affine_preimage(region, matrix, offset=None):
    """Return the exact preimage of ``region`` under ``x -> A*x+b``.

    Invertible square maps preserve canonical structure. Rectangular or singular
    maps are lowered by exact semialgebraic substitution.
    """
    ambient = (
        region.ambient_dimension
        if isinstance(region, HRepresentation)
        else region.ambient_dimension()
    )
    A = sp.ImmutableMatrix(matrix)
    if A.rows != ambient:
        raise ValueError("affine map row count must match the target ambient dimension")
    b = sp.ImmutableMatrix.zeros(A.rows, 1) if offset is None else sp.ImmutableMatrix(offset)
    if b.shape != (A.rows, 1):
        raise ValueError("affine offset must match the target ambient dimension")
    if _is_invertible_square(A):
        inv = A.inv()
        return affine_image(region, inv, tuple(-inv * b))

    source = collision_free_real_symbols("u", A.cols, region, A, b)
    mapping = tuple(A * sp.Matrix(source) + b)
    return region_preimage(region, mapping, source)


def _mapping_tuple(mapping):
    if isinstance(mapping, sp.Basic):
        return (sp.sympify(mapping),)
    return tuple(sp.sympify(expr) for expr in mapping)


def _mapping_variables(mapping, variables, expected_dimension: int):
    maps = _mapping_tuple(mapping)
    if variables is None:
        symbols = sorted(
            set().union(*(expr.free_symbols for expr in maps)), key=sp.default_sort_key
        )
        if len(symbols) != expected_dimension:
            raise ValueError(
                "variables are required when mapping symbols do not match the source dimension"
            )
        source = tuple(symbols)
    else:
        source = tuple(sp.Symbol(v, real=True) if isinstance(v, str) else v for v in variables)
    if len(source) != expected_dimension:
        raise ValueError("mapping variable count must match the source ambient dimension")
    return maps, source


def _affine_map_data(mapping, variables):
    maps = sp.Matrix(mapping)
    source = tuple(variables)
    A = maps.jacobian(source)
    if any(entry.free_symbols & set(source) for entry in A):
        return None
    b = sp.simplify(maps - A * sp.Matrix(source))
    if any(entry.free_symbols & set(source) for entry in b):
        return None
    return sp.ImmutableMatrix(A), sp.ImmutableMatrix(b)


def region_image(region, mapping, variables=None):
    """Return an exact symbolic image, preserving canonical affine structure when possible."""
    if isinstance(region, HRepresentation):
        source_dimension = region.ambient_dimension
    elif isinstance(region, StandardRegion):
        source_dimension = region.ambient_dimension()
    else:
        raise TypeError("region_image requires a canonical region or HRepresentation")
    maps, source = _mapping_variables(mapping, variables, source_dimension)
    affine = _affine_map_data(maps, source)
    if affine is not None:
        A, b = affine
        return affine_image(region, A, tuple(b))
    base = region.to_polytope() if isinstance(region, HRepresentation) else region
    return TransformedRegion(base, maps, source)


def region_preimage(target, mapping, variables=None, *, target_variables=None):
    """Return an exact symbolic preimage, preserving invertible affine structure when possible."""
    from .geometry_queries import semialgebraic_preimage
    from .symbolic_regions import SemialgebraicRegion, as_semialgebraic_region

    if isinstance(target, HRepresentation):
        target_dimension = target.ambient_dimension
    elif isinstance(target, StandardRegion):
        target_dimension = target.ambient_dimension()
    else:
        if not isinstance(target, SemialgebraicRegion):
            raise TypeError("region_preimage requires a canonical or semialgebraic target region")
        target_dimension = len(target.variables)

    maps = _mapping_tuple(mapping)
    if len(maps) != target_dimension:
        raise ValueError("mapping output dimension must match the target ambient dimension")
    if variables is None:
        symbols = sorted(
            set().union(*(expr.free_symbols for expr in maps)), key=sp.default_sort_key
        )
        if not symbols:
            raise ValueError("source variables are required for a constant preimage map")
        source = tuple(symbols)
    else:
        source = tuple(sp.Symbol(v, real=True) if isinstance(v, str) else v for v in variables)

    affine = _affine_map_data(maps, source)
    if (
        affine is not None
        and len(source) == target_dimension
        and _is_invertible_square(affine[0])
        and isinstance(target, (StandardRegion, HRepresentation))
    ):
        A, b = affine
        return affine_preimage(target, A, tuple(b))

    if isinstance(target, HRepresentation):
        target_vars = collision_free_real_symbols("y", target_dimension, target, maps, source)
        target_formula = target.as_formula(target_vars)
    elif isinstance(target, StandardRegion):
        target_vars = collision_free_real_symbols("y", target_dimension, target, maps, source)
        target_formula = as_semialgebraic_region(target, target_vars).quantifier_free_formula()
    else:
        target_vars = target.variables
        target_formula = target.quantifier_free_formula()
    if target_variables is not None:
        supplied = tuple(
            sp.Symbol(v, real=True) if isinstance(v, str) else v for v in target_variables
        )
        if len(supplied) != target_dimension:
            raise ValueError("target_variables must match the target ambient dimension")
        target_formula = target_formula.xreplace(dict(zip(target_vars, supplied, strict=True)))
        target_vars = supplied
    formula = semialgebraic_preimage(maps, target_formula, source, target_variables=target_vars)
    return SemialgebraicRegion(formula, source)


__all__ = ["affine_image", "affine_preimage", "region_image", "region_preimage"]
