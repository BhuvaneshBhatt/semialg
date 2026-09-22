"""Certified parametric descriptions of bounded semialgebraic geometry."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ._zero_testing import certified_equal
from .exact_arithmetic import compare_exact_reals
from .normalization import normalize_formula, normalize_problem_variables


@dataclass(frozen=True)
class ParametricChart:
    """One certified map from a bounded parameter domain into a region.

    ``condition`` is an additional condition on the parameters beyond the
    explicit bounds.  ``domain_dimension`` is present only when the parameter
    domain is certified to have that dimension; this lets dimension queries
    reuse a chart without running CAD.
    """

    parameters: tuple[sp.Symbol, ...]
    bounds: tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]
    condition: sp.Expr
    mapping: tuple[sp.Expr, ...]
    domain_dimension: int | None
    source: str

    @property
    def ambient_dimension(self) -> int:
        return len(self.mapping)

    def parameter_formula(self) -> sp.Expr:
        return normalize_formula(
            sp.And(
                self.condition,
                *(sp.And(param >= lo, param <= hi) for param, lo, hi in self.bounds),
            )
        )

    def jacobian_matrix(self) -> sp.ImmutableMatrix:
        """Return the exact ambient-by-parameter Jacobian of the chart."""

        if not self.parameters:
            return sp.ImmutableMatrix(sp.zeros(self.ambient_dimension, 0))
        return sp.ImmutableMatrix(sp.Matrix(self.mapping).jacobian(self.parameters))

    def gram_matrix(self) -> sp.ImmutableMatrix:
        """Return the pullback Euclidean metric ``J.T*J``."""

        jacobian = self.jacobian_matrix()
        return sp.ImmutableMatrix(jacobian.T * jacobian)

    def metric_factor(self) -> sp.Expr:
        """Return the intrinsic Hausdorff-measure density of this chart."""

        if not self.parameters:
            return sp.Integer(1)
        gram = self.gram_matrix()
        return sp.simplify(sp.sqrt(sp.det(gram)))

    def pullback(self, expression: object, variables: Sequence[sp.Symbol]) -> sp.Expr:
        """Pull an ambient expression back to this chart's parameter domain."""

        vars_ = tuple(map(sp.sympify, variables))
        if len(vars_) != self.ambient_dimension:
            raise ValueError("variable count does not match chart ambient dimension")
        return sp.simplify(sp.sympify(expression).subs(dict(zip(vars_, self.mapping, strict=True))))

    def intrinsic_integrand(self, expression: object, variables: Sequence[sp.Symbol]) -> sp.Expr:
        """Return the pullback integrand including the induced metric factor."""

        return sp.simplify(self.pullback(expression, variables) * self.metric_factor())

    def certified_image_dimension(self) -> int | None:
        """Return an exact image dimension when rank and domain data certify it."""

        if self.domain_dimension is None:
            return None
        if self.domain_dimension == 0:
            return 0
        if len(self.parameters) != self.domain_dimension:
            return None
        parameter_set = set(self.parameters)
        if any(expr.free_symbols - parameter_set for expr in self.mapping):
            return None
        try:
            rank = int(self.jacobian_matrix().rank())
        except (TypeError, ValueError, ArithmeticError, NotImplementedError):
            return None
        return min(rank, self.domain_dimension)


@dataclass(frozen=True)
class ParametricCover:
    """Finite certified cover of a region or a bounded region intersection."""

    charts: tuple[ParametricChart, ...]
    variables: tuple[sp.Symbol, ...]
    formula: sp.Expr
    exact: bool = True

    def certified_dimension(self) -> int | None:
        """Return the maximum certified chart-image dimension, if all are known."""

        if not self.charts:
            return -1 if self.formula in (False, sp.false) else None
        dims = tuple(chart.certified_image_dimension() for chart in self.charts)
        if any(dim is None for dim in dims):
            return None
        return max(int(dim) for dim in dims if dim is not None)


def _finite_bound(value: sp.Expr) -> bool:
    return value not in (-sp.oo, sp.oo) and value.is_finite is not False


def _bounds_full_dimensional(
    bounds: Sequence[tuple[sp.Symbol, sp.Expr, sp.Expr]],
) -> bool:
    for _param, lo, hi in bounds:
        if not (_finite_bound(lo) and _finite_bound(hi)):
            return False
        try:
            if compare_exact_reals(lo, hi) >= 0:
                return False
        except (TypeError, ValueError, NotImplementedError):
            return False
    return True


def _identity_chart(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    bounds: Sequence[tuple[sp.Symbol, sp.Expr, sp.Expr]],
    *,
    source: str,
) -> ParametricChart:
    box_formula = sp.And(*(sp.And(var >= lo, var <= hi) for var, lo, hi in bounds))
    condition = normalize_formula(sp.And(formula, box_formula))
    return ParametricChart(
        variables,
        tuple(bounds),
        condition,
        variables,
        None,
        source,
    )


def _normalize_clip_bounds(
    bounds: Sequence[Sequence[object]] | None,
    variables: tuple[sp.Symbol, ...],
) -> tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...] | None:
    if bounds is None:
        return None
    if len(bounds) != len(variables):
        raise ValueError("bounds must contain one (lower, upper) pair per variable")
    normalized: list[tuple[sp.Symbol, sp.Expr, sp.Expr]] = []
    for var, pair in zip(variables, bounds, strict=True):
        if len(pair) != 2:
            raise ValueError("each bound must be a (lower, upper) pair")
        lo, hi = map(sp.sympify, pair)
        try:
            if compare_exact_reals(lo, hi) > 0:
                raise ValueError(f"lower bound exceeds upper bound for {var}")
        except (TypeError, NotImplementedError):
            pass
        normalized.append((var, lo, hi))
    return tuple(normalized)


def _chart_from_parametric_region(region) -> ParametricChart | None:
    bounds = tuple((p, sp.sympify(lo), sp.sympify(hi)) for p, lo, hi in region.limits)
    if not all(_finite_bound(lo) and _finite_bound(hi) for _, lo, hi in bounds):
        return None
    condition = normalize_formula(region.assumptions)
    domain_dim = (
        len(region.parameters)
        if condition in (True, sp.true) and _bounds_full_dimensional(bounds)
        else None
    )
    return ParametricChart(
        tuple(region.parameters),
        bounds,
        condition,
        tuple(region.mapping),
        domain_dim,
        "parametric_region",
    )


def _simplex_chart(region, index: int) -> ParametricChart:
    vertices = tuple(region.vertices)
    dim = len(vertices) - 1
    if dim == 0:
        return ParametricChart((), (), sp.true, tuple(vertices[0]), 0, "simplex")
    params = tuple(sp.Dummy(f"u{index}_{i + 1}", real=True) for i in range(dim))
    origin = sp.Matrix(vertices[0])
    mapped = origin
    for param, vertex in zip(params, vertices[1:], strict=True):
        mapped += param * (sp.Matrix(vertex) - origin)
    # A cube bound plus the simplex condition gives a genuinely bounded chart
    # while avoiding dependent integration limits in the common chart model.
    bounds = tuple((param, sp.Integer(0), sp.Integer(1)) for param in params)
    condition = sp.Le(sum(params), 1)
    return ParametricChart(params, bounds, condition, tuple(mapped), dim, "simplex")


def _hyperspherical_coordinates(
    dimension: int, *, prefix: str
) -> tuple[
    tuple[sp.Symbol, ...], tuple[sp.Expr, ...], tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]
]:
    if dimension < 2:
        raise ValueError("hyperspherical coordinates require dimension >= 2")
    angles = tuple(sp.Dummy(f"{prefix}_theta_{i + 1}", real=True) for i in range(dimension - 1))
    coordinates: list[sp.Expr] = []
    sine_prefix = sp.Integer(1)
    for index in range(dimension):
        if index == dimension - 1:
            coordinate = sine_prefix
        else:
            coordinate = sine_prefix * sp.cos(angles[index])
            sine_prefix *= sp.sin(angles[index])
        coordinates.append(sp.simplify(coordinate))
    bounds = tuple(
        (angle, sp.Integer(0), 2 * sp.pi if index == dimension - 2 else sp.pi)
        for index, angle in enumerate(angles)
    )
    return angles, tuple(coordinates), bounds


def _radial_chart(
    center: Sequence[object],
    radius: object,
    *,
    surface: bool,
    affine: sp.MatrixBase | None = None,
    source: str,
) -> ParametricChart:
    center_vec = sp.ImmutableMatrix(tuple(map(sp.sympify, center)))
    dimension = len(center_vec)
    radius_expr = sp.sympify(radius)
    if dimension == 1:
        if surface:
            raise NotImplementedError(
                "one-dimensional sphere boundaries require a finite point cover"
            )
        radial = sp.Dummy("radial_u", real=True)
        return ParametricChart(
            (radial,),
            ((radial, -radius_expr, radius_expr),),
            sp.true,
            (center_vec[0] + radial,),
            1,
            source,
        )
    angles, omega_values, angular_bounds = _hyperspherical_coordinates(dimension, prefix=source)
    omega = sp.ImmutableMatrix(omega_values)
    linear = sp.eye(dimension) if affine is None else sp.ImmutableMatrix(affine)
    if linear.shape != (dimension, dimension):
        raise ValueError("radial affine map must match the ambient dimension")
    if surface:
        parameters = angles
        bounds = angular_bounds
        mapped = center_vec + radius_expr * linear * omega
        domain_dimension = dimension - 1
    else:
        radial = sp.Dummy(f"{source}_r", nonnegative=True, real=True)
        parameters = (radial, *angles)
        bounds = ((radial, sp.Integer(0), radius_expr), *angular_bounds)
        mapped = center_vec + radial * linear * omega
        domain_dimension = dimension
    return ParametricChart(
        parameters,
        tuple(bounds),
        sp.true,
        tuple(map(sp.simplify, mapped)),
        domain_dimension,
        source,
    )


def _ellipsoid_factor(shape: sp.MatrixBase) -> sp.ImmutableMatrix:
    matrix = sp.ImmutableMatrix(shape)
    try:
        return sp.ImmutableMatrix(matrix.cholesky())
    except (ValueError, TypeError, NotImplementedError) as exc:
        raise NotImplementedError(
            "could not construct an exact positive-definite ellipsoid factor"
        ) from exc


def _standard_region_charts(region) -> tuple[ParametricChart, ...] | None:
    from .standard_regions import (
        Ball,
        Box,
        Ellipsoid,
        EllipsoidBoundary,
        FinitePointSet,
        Interval,
        Parallelepiped,
        Parallelogram,
        ParametricRegion,
        Polygon,
        Simplex,
        Sphere,
        TetrahedralComplex,
    )

    if isinstance(region, ParametricRegion):
        chart = _chart_from_parametric_region(region)
        return None if chart is None else (chart,)
    if isinstance(region, Sphere):
        if region.ambient_dimension() == 1:
            left = tuple(sp.simplify(c - region.radius) for c in region.center)
            right = tuple(sp.simplify(c + region.radius) for c in region.center)
            return (
                ParametricChart((), (), sp.true, left, 0, "sphere_point"),
                ParametricChart((), (), sp.true, right, 0, "sphere_point"),
            )
        return (_radial_chart(region.center, region.radius, surface=True, source="sphere"),)
    if isinstance(region, Ball):
        return (_radial_chart(region.center, region.radius, surface=False, source="ball"),)
    if isinstance(region, EllipsoidBoundary):
        factor = _ellipsoid_factor(region.shape_matrix)
        return (
            _radial_chart(
                region.center, 1, surface=True, affine=factor, source="ellipsoid_boundary"
            ),
        )
    if isinstance(region, Ellipsoid):
        factor = _ellipsoid_factor(region.shape_matrix)
        return (_radial_chart(region.center, 1, surface=False, affine=factor, source="ellipsoid"),)
    if isinstance(region, FinitePointSet):
        return tuple(
            ParametricChart((), (), sp.true, tuple(point), 0, "point") for point in region.points
        )
    if isinstance(region, Interval):
        if not (_finite_bound(region.lower) and _finite_bound(region.upper)):
            return None
        try:
            cmp = compare_exact_reals(region.lower, region.upper)
        except (TypeError, ValueError, NotImplementedError):
            cmp = 0 if certified_equal(region.upper, region.lower) is True else None
        if cmp == 0:
            if not (region.lower_closed and region.upper_closed):
                return ()
            return (ParametricChart((), (), sp.true, (region.lower,), 0, "interval"),)
        param = sp.Dummy("u", real=True)
        condition = sp.true
        if not region.lower_closed:
            condition = sp.And(condition, param > region.lower)
        if not region.upper_closed:
            condition = sp.And(condition, param < region.upper)
        return (
            ParametricChart(
                (param,),
                ((param, region.lower, region.upper),),
                condition,
                (param,),
                1 if cmp == -1 else None,
                "interval",
            ),
        )
    if isinstance(region, Box):
        params: list[sp.Symbol] = []
        bounds: list[tuple[sp.Symbol, sp.Expr, sp.Expr]] = []
        mapping: list[sp.Expr] = []
        unknown_dimension = False
        for i, (raw_lo, raw_hi) in enumerate(region.bounds):
            lo, hi = sp.sympify(raw_lo), sp.sympify(raw_hi)
            if not (_finite_bound(lo) and _finite_bound(hi)):
                return None
            try:
                cmp = compare_exact_reals(lo, hi)
            except (TypeError, ValueError, NotImplementedError):
                cmp = 0 if certified_equal(hi, lo) is True else None
            if cmp == 0:
                mapping.append(lo)
                continue
            param = sp.Dummy(f"u{i + 1}", real=True)
            params.append(param)
            bounds.append((param, lo, hi))
            mapping.append(param)
            if cmp is None:
                unknown_dimension = True
        return (
            ParametricChart(
                tuple(params),
                tuple(bounds),
                sp.true,
                tuple(mapping),
                None if unknown_dimension else len(params),
                "box",
            ),
        )
    if isinstance(region, Simplex):
        return (_simplex_chart(region, 0),)
    if isinstance(region, Polygon):
        return tuple(
            _simplex_chart(triangle, i) for i, triangle in enumerate(region.triangulation())
        )
    if isinstance(region, TetrahedralComplex):
        return tuple(_simplex_chart(tet, i) for i, tet in enumerate(region.tetrahedra))
    if isinstance(region, (Parallelogram, Parallelepiped)):
        params = tuple(sp.Dummy(f"u{i + 1}", real=True) for i in range(len(region.vectors)))
        mapped = sp.Matrix(region.origin)
        for param, vector in zip(params, region.vectors, strict=True):
            mapped += param * sp.Matrix(vector)
        bounds = tuple((param, sp.Integer(0), sp.Integer(1)) for param in params)
        return (
            ParametricChart(
                params,
                bounds,
                sp.true,
                tuple(mapped),
                region.dimension(),
                "parallelepiped",
            ),
        )
    return None


def bounded_parametric_cover(
    region: object,
    variables: Sequence[sp.Symbol | str] | None = None,
    bounds: Sequence[Sequence[object]] | None = None,
) -> ParametricCover:
    """Return a certified finite bounded parametric cover when one is structural.

    For recognized bounded ``StandardRegion`` objects, the cover uses their
    natural low-dimensional parametrizations.  When explicit clipping bounds
    are supplied, an identity chart over the bounded box represents the exact
    intersection without invoking CAD.  Unsupported unbounded inputs raise
    ``NotImplementedError`` rather than silently requesting decomposition.
    """

    from .standard_regions import StandardRegion
    from .symbolic_regions import as_semialgebraic_region

    if isinstance(region, StandardRegion):
        ambient = region.ambient_dimension()
        vars_ = tuple(
            normalize_problem_variables(
                variables
                if variables is not None
                else tuple(sp.Symbol(f"x{i + 1}", real=True) for i in range(ambient)),
                sp.true,
            )
        )
        if len(vars_) != ambient:
            raise ValueError("variable count does not match region ambient dimension")
        symbolic = as_semialgebraic_region(region, vars_)
        clip_bounds = _normalize_clip_bounds(bounds, vars_)
        if clip_bounds is not None:
            chart = _identity_chart(
                symbolic.formula, vars_, clip_bounds, source="bounded_intersection"
            )
            return ParametricCover((chart,), vars_, chart.parameter_formula(), True)
        charts = _standard_region_charts(region)
        if charts is None:
            raise NotImplementedError(
                "no bounded structural parametrization is available for this region"
            )
        return ParametricCover(tuple(charts), vars_, symbolic.formula, True)

    symbolic = as_semialgebraic_region(region, variables)
    vars_ = tuple(symbolic.variables)
    clip_bounds = _normalize_clip_bounds(bounds, vars_)
    if clip_bounds is None:
        raise NotImplementedError(
            "formula regions require explicit finite bounds for a bounded parametric cover"
        )
    if not all(_finite_bound(lo) and _finite_bound(hi) for _, lo, hi in clip_bounds):
        raise ValueError("bounded parametric covers require finite clipping bounds")
    formula = symbolic.formula
    chart = _identity_chart(formula, vars_, clip_bounds, source="bounded_intersection")
    return ParametricCover((chart,), vars_, chart.parameter_formula(), True)


def intrinsic_parametric_cover(
    region: object,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    dimension: int | None = None,
    require_verified: bool = True,
) -> ParametricCover:
    """Return exact charts for intrinsic geometry, using CAD when necessary.

    Canonical ``StandardRegion`` objects use structural parametrizations first.
    Formula regions are decomposed into certified regular CAD strata and those
    graph charts are normalized to :class:`ParametricChart`.  Singular target-
    dimensional strata are never silently discarded.
    """

    from .standard_regions import StandardRegion

    if isinstance(region, StandardRegion):
        cover = bounded_parametric_cover(region, variables)
        expected = region.dimension() if dimension is None else int(dimension)
        if expected != region.dimension():
            raise ValueError(
                "requested intrinsic dimension does not match canonical region dimension"
            )
        for chart in cover.charts:
            measured = chart.certified_image_dimension()
            if measured is not None and measured != expected:
                raise ValueError(
                    f"structural chart dimension {measured} does not match region dimension {expected}"
                )
        return cover

    from .cad_algorithms.cells import (
        extract_cylindrical_solution,
        intrinsic_cell_integral,
        stratify_intrinsic_solution,
    )
    from .regions.operations import region_dimension
    from .symbolic_regions import as_semialgebraic_region

    symbolic = as_semialgebraic_region(region, variables)
    vars_ = tuple(symbolic.variables)
    formula = symbolic.quantifier_free_formula()
    target = region_dimension(formula, vars_) if dimension is None else int(dimension)
    if target < 0 or target > len(vars_):
        raise ValueError("intrinsic dimension must be between 0 and the ambient dimension")
    solution = extract_cylindrical_solution(formula, vars_, selected_only=True)
    stratification = stratify_intrinsic_solution(
        solution, dimension=target, require_verified=require_verified
    )
    singular = tuple(stratification.singular_strata)
    if singular:
        indices = tuple(tuple(stratum.cell_index) for stratum in singular)
        raise ValueError(f"intrinsic CAD cover contains uncertified singular strata: {indices}")

    charts: list[ParametricChart] = []
    for stratum in stratification.regular_strata:
        piece = intrinsic_cell_integral(
            stratum.cell, sp.Integer(1), evaluate=False, require_verified=require_verified
        )
        mapping_dict = dict(piece.chart_mapping)
        mapping = tuple(sp.simplify(mapping_dict[var]) for var in vars_)
        # CAD integral limits are stored in integration order; chart bounds are
        # coordinate-order metadata, so restore parameter order here.
        limit_map = {var: (sp.sympify(lo), sp.sympify(hi)) for var, lo, hi in piece.limits}
        params = tuple(piece.chart_variables)
        bounds = tuple((var, *limit_map[var]) for var in params)
        charts.append(
            ParametricChart(
                params,
                bounds,
                sp.true,
                mapping,
                target,
                "cad_intrinsic",
            )
        )
    return ParametricCover(tuple(charts), vars_, formula, True)


def certify_parametric_dimension(
    region: object,
    variables: Sequence[sp.Symbol | str] | None = None,
) -> int | None:
    """Certify region dimension from structural charts without invoking CAD."""

    try:
        cover = bounded_parametric_cover(region, variables)
    except (TypeError, ValueError, NotImplementedError):
        return None
    return cover.certified_dimension()


__all__ = [
    "ParametricChart",
    "ParametricCover",
    "bounded_parametric_cover",
    "intrinsic_parametric_cover",
    "certify_parametric_dimension",
]
