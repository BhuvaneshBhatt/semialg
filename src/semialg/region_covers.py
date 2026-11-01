"""Structural covers and exact outer boxes for geometric optimization."""

from dataclasses import fields, is_dataclass, replace

import sympy as sp

from .parametric_geometry import ParametricCover, bounded_parametric_cover
from .standard_regions import (
    Ball,
    BooleanRegion,
    Box,
    Ellipsoid,
    EllipsoidBoundary,
    FinitePointSet,
    Interval,
    ParametricRegion,
    Simplex,
    Sphere,
    TransformedRegion,
)


def _inexact_data(value):
    if isinstance(value, (sp.Basic, sp.MatrixBase)):
        return value.has(sp.Float)
    if is_dataclass(value):
        return any(_inexact_data(getattr(value, field.name)) for field in fields(value))
    if isinstance(value, (tuple, list)):
        return any(_inexact_data(item) for item in value)
    return isinstance(value, float)


def optimization_cover(region, variables):
    """Return a covering parameterization, retaining boundaries and point strata.

    Structural charts cover the entire region; they are not the regular,
    top-dimensional charts used for integration. Unions may overlap. Affine or
    nonlinear transformed regions compose the base cover without inversion.
    """
    if _inexact_data(region):
        raise NotImplementedError("structural optimization geometry requires exact coefficients")
    coordinates = tuple(variables)
    if isinstance(region, BooleanRegion) and region.op == "union":
        covers = tuple(optimization_cover(part, coordinates) for part in region.regions)
        return ParametricCover(
            tuple(chart for cover in covers for chart in cover.charts),
            coordinates,
            region.as_formula(coordinates),
            True,
        )
    if isinstance(region, TransformedRegion):
        cover = optimization_cover(region.base, region.base_variables)
        charts = tuple(
            replace(
                chart,
                mapping=tuple(
                    expression.xreplace(
                        dict(zip(region.base_variables, chart.mapping, strict=True))
                    )
                    for expression in region.mapping
                ),
                source="transformed:" + chart.source,
            )
            for chart in cover.charts
        )
        return ParametricCover(charts, coordinates, region.as_formula(coordinates), True)
    cover = bounded_parametric_cover(region, coordinates)
    return replace(
        cover,
        charts=tuple(
            replace(chart, condition=sp.And(chart.condition, *region.conditions))
            for chart in cover.charts
        ),
    )


def _interval(expression, bounds):
    """Enclose an exact expression using dependency-resolved chart bounds.

    Interval operations may overestimate ranges. Even powers include zero if
    the input interval straddles it; unsupported maps are not guessed.
    """
    if expression in bounds:
        return bounds[expression]
    if not expression.free_symbols.intersection(bounds):
        return expression, expression
    if expression.is_Add:
        parts = [_interval(arg, bounds) for arg in expression.args]
        return sum(p[0] for p in parts), sum(p[1] for p in parts)
    if expression.is_Mul:
        lo = hi = sp.Integer(1)
        for arg in expression.args:
            a, b = _interval(arg, bounds)
            products = (lo * a, lo * b, hi * a, hi * b)
            lo, hi = sp.Min(*products), sp.Max(*products)
        return lo, hi
    if expression.is_Pow and expression.exp.is_Integer and expression.exp >= 0:
        lo, hi = _interval(expression.base, bounds)
        power = expression.exp
        if power % 2:
            return lo**power, hi**power
        return sp.Piecewise(
            (0, sp.And(lo <= 0, hi >= 0)), (sp.Min(lo**power, hi**power), True)
        ), sp.Max(lo**power, hi**power)
    raise NotImplementedError("no structural interval rule for this coordinate map")


def _parameter_bounds(limits):
    """Resolve nested chart limits in dependency order, independent of tuple order."""
    variables = {variable for variable, _, _ in limits}
    pending = list(limits)
    bounds = {}
    while pending:
        ready = [
            (v, lo, hi)
            for v, lo, hi in pending
            if ((lo.free_symbols | hi.free_symbols) & variables) <= bounds.keys()
        ]
        if not ready:
            raise NotImplementedError("cyclic chart limits have no structural outer box")
        for variable, lower, upper in ready:
            bounds[variable] = (_interval(lower, bounds)[0], _interval(upper, bounds)[1])
            pending.remove((variable, lower, upper))
    return bounds


def bounding_cover(region, variables):
    """Return an exact outer coordinate box without calling an optimizer or CAD.

    Bounds need not be tight. Missing structural rules raise NotImplementedError.
    Symbolic endpoints retain parameters and are not numerical search bounds
    until specialized. Membership in the original region remains necessary.
    """
    if _inexact_data(region):
        raise NotImplementedError("structural optimization geometry requires exact coefficients")
    coordinates = tuple(variables)
    if len(coordinates) != region.ambient_dimension():
        raise ValueError("variable count does not match region ambient dimension")
    if isinstance(region, (Ball, Sphere)):
        pairs = tuple((c - region.radius, c + region.radius) for c in region.center)
    elif isinstance(region, (Ellipsoid, EllipsoidBoundary)):
        pairs = tuple(
            (c - sp.sqrt(region.shape_matrix[i, i]), c + sp.sqrt(region.shape_matrix[i, i]))
            for i, c in enumerate(region.center)
        )
    elif isinstance(region, Box):
        pairs = region.bounds
    elif isinstance(region, Interval):
        pairs = ((region.lower, region.upper),)
    elif isinstance(region, (Simplex, FinitePointSet)):
        points = region.vertices if isinstance(region, Simplex) else region.points
        if not points:
            raise NotImplementedError("an empty point set has no coordinate extrema")
        pairs = tuple(
            (sp.Min(*(p[i] for p in points)), sp.Max(*(p[i] for p in points)))
            for i in range(len(coordinates))
        )
    elif isinstance(region, BooleanRegion):
        if not region.regions or region.op in {"complement", "symmetric_difference"}:
            raise NotImplementedError("no structural outer box for this Boolean region")
        boxes = tuple(bounding_cover(part, coordinates) for part in region.regions)
        if region.op == "union":
            pairs = tuple(
                (sp.Min(*(box[v][0] for box in boxes)), sp.Max(*(box[v][1] for box in boxes)))
                for v in coordinates
            )
        else:
            # Either side covers an intersection; only the left covers a difference.
            return boxes[0]
    elif isinstance(region, TransformedRegion):
        bounds = bounding_cover(region.base, region.base_variables)
        pairs = tuple(_interval(expr, bounds) for expr in region.mapping)
    elif isinstance(region, ParametricRegion):
        bounds = _parameter_bounds(region.limits)
        pairs = tuple(_interval(expr, bounds) for expr in region.mapping)
    else:
        cover = optimization_cover(region, coordinates)
        boxes = []
        for chart in cover.charts:
            bounds = _parameter_bounds(chart.bounds)
            boxes.append(tuple(_interval(expr, bounds) for expr in chart.mapping))
        if not boxes:
            raise NotImplementedError("no nonempty structural outer box")
        pairs = tuple(
            (sp.Min(*(box[i][0] for box in boxes)), sp.Max(*(box[i][1] for box in boxes)))
            for i in range(len(coordinates))
        )
    return dict(zip(coordinates, pairs, strict=True))
