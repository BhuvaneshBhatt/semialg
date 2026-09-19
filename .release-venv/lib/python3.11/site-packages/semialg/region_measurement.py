from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from .measure import MeasureResult
from .moments import RegionCentroidResult


def _variables_for_region(region) -> tuple[sp.Symbol, ...]:
    return tuple(sp.Symbol(f"x{i + 1}", real=True) for i in range(region.ambient_dimension()))


def _target_dimension(region, measure_dimension: object) -> int:
    if measure_dimension in (None, "intrinsic"):
        return region.dimension()
    if measure_dimension == "ambient":
        return region.ambient_dimension()
    try:
        target = int(measure_dimension)
    except (TypeError, ValueError) as exc:
        raise ValueError("measure_dimension must be 'intrinsic', 'ambient', or an integer") from exc
    if target < 0 or target > region.ambient_dimension():
        raise ValueError("measure dimension is outside the ambient dimension")
    return target


def _ellipsoid_volume(region) -> sp.Expr:
    n = region.ambient_dimension()
    unit = sp.pi ** sp.Rational(n, 2) / sp.gamma(sp.Rational(n, 2) + 1)
    return sp.simplify(sp.sqrt(sp.det(region.shape_matrix)) * unit)


def _chart_integral(chart, expression: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr:
    integrand = chart.intrinsic_integrand(expression, variables)
    if chart.condition not in (True, sp.true):
        from .region_integrate import integrate_over_region

        return sp.sympify(
            integrate_over_region(
                integrand,
                chart.parameter_formula(),
                chart.parameters,
                measure_dimension="ambient",
            )
        )
    value = integrand
    for parameter, lower, upper in reversed(chart.bounds):
        value = sp.integrate(value, (parameter, lower, upper))
    return sp.simplify(value)


def region_measure(
    region: object,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    measure_dimension: object = "intrinsic",
    return_result: bool = False,
) -> sp.Expr | MeasureResult:
    """Return exact Euclidean/Hausdorff measure of a region when supported.

    Canonical geometry defaults to intrinsic measure. Formula regions delegate
    to :func:`semialgebraic_measure` and therefore retain its ambient-measure
    behavior unless another dimension is requested explicitly.
    """

    from .standard_regions import Ellipsoid, StandardRegion

    if not isinstance(region, StandardRegion):
        from .measure import semialgebraic_measure

        if variables is None:
            raise ValueError("variables are required for formula-region measure")
        return semialgebraic_measure(
            region,
            variables,
            measure_dimension=measure_dimension,
            return_result=return_result,
        )

    dimension = region.dimension()
    ambient = region.ambient_dimension()
    target = _target_dimension(region, measure_dimension)
    vars_ = _variables_for_region(region)

    if dimension < 0:
        value = sp.Integer(0)
        method = "empty_region"
    elif target > dimension:
        value = sp.Integer(0)
        method = "dimension_comparison"
    elif target < dimension:
        value = sp.oo
        method = "dimension_comparison"
    elif isinstance(region, Ellipsoid):
        value = _ellipsoid_volume(region)
        method = "ellipsoid_determinant"
    else:
        from .standard_region_integrate import integrate_over_standard_region

        try:
            value = sp.simplify(integrate_over_standard_region(1, region, vars_))
            method = "standard_region_integration"
        except NotImplementedError:
            try:
                cover = region.intrinsic_parametric_cover(vars_)
                value = sp.simplify(
                    sp.Add(
                        *(_chart_integral(chart, sp.Integer(1), vars_) for chart in cover.charts)
                    )
                )
                method = "intrinsic_parametric_cover"
            except NotImplementedError:
                from .region_integrate import integrate_over_region

                value = sp.sympify(
                    integrate_over_region(
                        1,
                        region.as_formula(vars_, eliminate=True),
                        vars_,
                        measure_dimension="intrinsic",
                    )
                )
                method = "semialgebraic_integration"

    result = MeasureResult(
        sp.sympify(value),
        vars_,
        method,
        {
            "ambient_dimension": ambient,
            "intrinsic_dimension": dimension,
            "measure_dimension": target,
        },
    )
    return result if return_result else result.value


def geometry_centroid(
    region: object,
    *,
    measure_dimension: object = "intrinsic",
    return_result: bool = False,
):
    """Return the centroid of a canonical geometry under uniform intrinsic measure."""

    from .standard_regions import (
        BallRegion,
        Ellipsoid,
        EllipsoidBoundary,
        ParallelepipedRegion,
        ParallelogramRegion,
        PointRegion,
        SimplexRegion,
        SphereRegion,
        SphericalShellRegion,
        StandardRegion,
    )

    if not isinstance(region, StandardRegion):
        raise TypeError("geometry_centroid expects a StandardRegion")
    dimension = region.dimension()
    target = _target_dimension(region, measure_dimension)
    if target != dimension:
        raise ValueError("centroid requires the region's intrinsic measure dimension")
    vars_ = _variables_for_region(region)
    measure = sp.sympify(region_measure(region, measure_dimension=target))
    if measure == 0 or measure in (sp.oo, -sp.oo, sp.zoo):
        raise ValueError("centroid requires finite nonzero measure")

    if isinstance(region, PointRegion):
        coords = tuple(
            sp.simplify(sum(point[i] for point in region.points) / len(region.points))
            for i in range(region.ambient_dimension())
        )
        method = "point_average"
    elif isinstance(
        region, (BallRegion, SphereRegion, SphericalShellRegion, Ellipsoid, EllipsoidBoundary)
    ):
        coords = tuple(region.center)
        method = "central_symmetry"
    elif isinstance(region, SimplexRegion):
        coords = tuple(
            sp.simplify(sum(vertex[i] for vertex in region.vertices) / len(region.vertices))
            for i in range(region.ambient_dimension())
        )
        method = "simplex_barycenter"
    elif isinstance(region, (ParallelogramRegion, ParallelepipedRegion)):
        coords = tuple(
            sp.simplify(
                region.origin[i] + sp.Rational(1, 2) * sum(vector[i] for vector in region.vectors)
            )
            for i in range(region.ambient_dimension())
        )
        method = "central_symmetry"
    else:
        from .standard_region_integrate import integrate_over_standard_region

        try:
            coords = tuple(
                sp.simplify(integrate_over_standard_region(var, region, vars_) / measure)
                for var in vars_
            )
            method = "standard_region_integration"
        except NotImplementedError:
            try:
                cover = region.intrinsic_parametric_cover(vars_)
                coords = tuple(
                    sp.simplify(
                        sp.Add(*(_chart_integral(chart, var, vars_) for chart in cover.charts))
                        / measure
                    )
                    for var in vars_
                )
                method = "intrinsic_parametric_cover"
            except NotImplementedError:
                from .region_integrate import integrate_over_region

                formula = region.as_formula(vars_, eliminate=True)
                coords = tuple(
                    sp.simplify(
                        sp.sympify(
                            integrate_over_region(
                                var, formula, vars_, measure_dimension="intrinsic"
                            )
                        )
                        / measure
                    )
                    for var in vars_
                )
                method = "semialgebraic_integration"

    mapping: Mapping[sp.Symbol, sp.Expr] = dict(zip(vars_, coords, strict=True))
    result = RegionCentroidResult(
        centroid=mapping,
        measure=measure,
        condition=region.as_formula(vars_, eliminate=True),
        variables=vars_,
        method=method,
        diagnostics={"measure_dimension": target},
        exact=not any(sp.sympify(value).has(sp.Float) for value in coords),
    )
    return result if return_result else coords


__all__ = ["geometry_centroid", "region_measure"]
