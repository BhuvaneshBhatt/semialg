from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp


@dataclass(frozen=True)
class SolutionPlotData:
    """Lightweight geometric data extracted from a semialgebraic solution.

    The object intentionally avoids depending on plotting libraries. Points,
    segments, and polygons are enough for tests, downstream adapters, and simple
    Matplotlib rendering through :func:`plot_solution`.
    """

    variables: tuple[sp.Symbol, ...]
    points: tuple[tuple[sp.Expr, ...], ...] = ()
    segments: tuple[tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]], ...] = ()
    polygons: tuple[tuple[tuple[sp.Expr, ...], ...], ...] = ()
    samples: tuple[Mapping[sp.Symbol, sp.Expr], ...] = ()
    source: str = "unknown"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    @property
    def dimension(self) -> int:
        if self.polygons:
            return 2
        if self.segments:
            return 1
        if self.points:
            return 0
        return -1

    def to_dict(self) -> dict[str, object]:
        return {
            "variables": tuple(sp.sstr(v) for v in self.variables),
            "points": tuple(tuple(sp.sstr(c) for c in point) for point in self.points),
            "segments": tuple(
                (tuple(sp.sstr(c) for c in start), tuple(sp.sstr(c) for c in end))
                for start, end in self.segments
            ),
            "polygons": tuple(
                tuple(tuple(sp.sstr(c) for c in point) for point in polygon)
                for polygon in self.polygons
            ),
            "sample_count": len(self.samples),
            "source": self.source,
            "diagnostics": dict(self.diagnostics),
        }


def _numeric_pair(point: Sequence[sp.Expr]) -> tuple[float, float]:
    if len(point) != 2:
        raise ValueError("expected a 2D point")
    return (float(sp.N(point[0])), float(sp.N(point[1])))


def _interval_endpoint(value: sp.Expr, default: sp.Expr) -> sp.Expr:
    return default if value in {-sp.oo, sp.oo} else value


def discretize_solution(
    solution,
    *,
    bounds: Sequence[tuple[sp.Expr, sp.Expr]] | None = None,
    samples_per_curve: int = 33,
) -> SolutionPlotData:
    """Return a small plotting/discretization representation for a solution.

    Supported exact inputs are one-dimensional interval components and 2D
    vertical cells. For curved 2D bounds, boundary curves are sampled at a
    modest deterministic grid. The output is intended for visualization and
    debugging, not as a certified mesh.
    """

    variables = tuple(getattr(solution, "variables", ()))
    bounds = tuple(bounds or ())
    points: list[tuple[sp.Expr, ...]] = []
    segments: list[tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]]] = []
    polygons: list[tuple[tuple[sp.Expr, ...], ...]] = []

    components = tuple(getattr(solution, "components", ()) or ())
    cells = tuple(getattr(solution, "cells", ()) or ())
    samples = tuple(getattr(solution, "samples", ()) or ())

    if len(variables) == 1 and components:
        default_lo, default_hi = (-sp.Integer(10), sp.Integer(10))
        if bounds:
            default_lo, default_hi = bounds[0]
        for component in components:
            lower = _interval_endpoint(getattr(component, "lower", -sp.oo), default_lo)
            upper = _interval_endpoint(getattr(component, "upper", sp.oo), default_hi)
            if getattr(component, "is_point", False):
                points.append((lower,))
            else:
                segments.append(((lower,), (upper,)))
        return SolutionPlotData(
            variables=variables,
            points=tuple(points),
            segments=tuple(segments),
            samples=samples,
            source="interval-components",
            diagnostics={"component_count": len(components)},
        )

    if len(variables) == 2 and cells:
        x, y = variables
        for cell in cells:
            x_interval = getattr(cell, "x_interval", None)
            y_bounds = tuple(getattr(cell, "y_bounds", ()) or ())
            if x_interval is None or not y_bounds:
                sample_fn = getattr(cell, "sample_point", None)
                if sample_fn is not None:
                    sample = sample_fn()
                    if isinstance(sample, Mapping) and x in sample and y in sample:
                        points.append((sample[x], sample[y]))
                continue
            xlo, xhi = x_interval
            if xlo in {-sp.oo, sp.oo} or xhi in {-sp.oo, sp.oo}:
                if bounds and len(bounds) >= 1:
                    bxlo, bxhi = bounds[0]
                    xlo = bxlo if xlo == -sp.oo else xlo
                    xhi = bxhi if xhi == sp.oo else xhi
                else:
                    # Unbounded 2D cells need explicit plotting bounds.
                    continue
            if xlo == xhi:
                for lower, upper in y_bounds:
                    lower_v = sp.simplify(lower.subs(x, xlo)) if hasattr(lower, "subs") else lower
                    upper_v = sp.simplify(upper.subs(x, xlo)) if hasattr(upper, "subs") else upper
                    if lower_v == upper_v:
                        points.append((xlo, lower_v))
                    else:
                        segments.append(((xlo, lower_v), (xlo, upper_v)))
                continue
            n = max(2, int(samples_per_curve))
            xs = [sp.simplify(xlo + (xhi - xlo) * sp.Rational(i, n - 1)) for i in range(n)]
            for lower, upper in y_bounds:
                lower_curve = [
                    (xv, sp.simplify(lower.subs(x, xv)) if hasattr(lower, "subs") else lower)
                    for xv in xs
                ]
                upper_curve = [
                    (xv, sp.simplify(upper.subs(x, xv)) if hasattr(upper, "subs") else upper)
                    for xv in xs
                ]
                if all(sp.simplify(upper_curve[i][1] - lower_curve[i][1]) == 0 for i in range(n)):
                    segments.extend((lower_curve[i], lower_curve[i + 1]) for i in range(n - 1))
                else:
                    polygons.append(tuple(lower_curve + list(reversed(upper_curve))))
        return SolutionPlotData(
            variables=variables,
            points=tuple(points),
            segments=tuple(segments),
            polygons=tuple(polygons),
            samples=samples,
            source="vertical-cells-2d",
            diagnostics={"cell_count": len(cells), "samples_per_curve": samples_per_curve},
        )

    # Fallback: return stored samples as points if their coordinates match the
    # requested variables.
    for sample in samples:
        if all(var in sample for var in variables):
            points.append(tuple(sample[var] for var in variables))
    return SolutionPlotData(
        variables=variables,
        points=tuple(points),
        samples=samples,
        source="samples",
        diagnostics={"sample_count": len(samples)},
    )


def _has_drawable_geometry(data: SolutionPlotData) -> bool:
    return bool(data.points or data.segments or data.polygons)


def _formula_from_solution(solution) -> sp.Expr | None:
    formula = getattr(solution, "formula", None)
    return formula if formula is not None else None


def _truth_mask_from_formula(formula: sp.Expr, variables: Sequence[sp.Symbol], xs, ys):
    """Evaluate a 2D Boolean formula on a NumPy mesh."""

    import numpy as np

    x_var, y_var = variables
    xx, yy = np.meshgrid(xs, ys)
    try:
        func = sp.lambdify((x_var, y_var), formula, modules="numpy")
        values = func(xx, yy)
        if np.isscalar(values):
            return np.full(xx.shape, bool(values), dtype=bool)
        return np.asarray(values, dtype=bool)
    except Exception:
        mask = np.zeros(xx.shape, dtype=bool)
        for row in range(xx.shape[0]):
            for col in range(xx.shape[1]):
                try:
                    value = formula.subs({x_var: float(xx[row, col]), y_var: float(yy[row, col])})
                    mask[row, col] = bool(value)
                except Exception:
                    mask[row, col] = False
        return mask


def _plot_formula_region_fallback(
    ax,
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    bounds,
    *,
    resolution: int,
    alpha: float,
    **plot_kwargs,
) -> bool:
    """Draw a sampled filled 2D region for a raw Boolean formula."""

    if len(variables) != 2 or not bounds or len(bounds) < 2:
        return False
    import numpy as np

    xlo, xhi = float(sp.N(bounds[0][0])), float(sp.N(bounds[0][1]))
    ylo, yhi = float(sp.N(bounds[1][0])), float(sp.N(bounds[1][1]))
    if not all(np.isfinite([xlo, xhi, ylo, yhi])) or xlo >= xhi or ylo >= yhi:
        return False
    n = max(16, int(resolution))
    xs = np.linspace(xlo, xhi, n)
    ys = np.linspace(ylo, yhi, n)
    mask = _truth_mask_from_formula(formula, variables, xs, ys)
    if not np.any(mask):
        return False
    xx, yy = np.meshgrid(xs, ys)
    ax.contourf(xx, yy, mask.astype(float), levels=[0.5, 1.5], alpha=alpha, **plot_kwargs)
    return True


def plot_solution(
    solution,
    *,
    bounds: Sequence[tuple[sp.Expr, sp.Expr]] | None = None,
    samples_per_curve: int = 33,
    raster_resolution: int = 300,
    ax=None,
    show: bool = False,
    **plot_kwargs,
):
    """Plot a 1D/2D solution using Matplotlib when available.

    The function returns the Matplotlib axes object. It is a convenience hook
    over :func:`discretize_solution`; callers that need certified geometry
    should use the returned discretization data instead.
    """

    data = discretize_solution(solution, bounds=bounds, samples_per_curve=samples_per_curve)
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency path
        raise ImportError("plot_solution requires matplotlib") from exc

    if ax is None:
        _, ax = plt.subplots()
    if len(data.variables) == 1:
        for point in data.points:
            ax.plot([float(sp.N(point[0]))], [0.0], marker="o", linestyle="None", **plot_kwargs)
        for start, end in data.segments:
            ax.plot([float(sp.N(start[0])), float(sp.N(end[0]))], [0.0, 0.0], **plot_kwargs)
        ax.set_xlabel(sp.sstr(data.variables[0]))
        ax.set_yticks([])
    elif len(data.variables) == 2:
        alpha = plot_kwargs.pop("alpha", 0.25)
        if not _has_drawable_geometry(data):
            formula = _formula_from_solution(solution)
            if formula is not None:
                _plot_formula_region_fallback(
                    ax,
                    formula,
                    data.variables,
                    bounds,
                    resolution=raster_resolution,
                    alpha=alpha,
                    **plot_kwargs,
                )
        for polygon in data.polygons:
            xs = [float(sp.N(point[0])) for point in polygon]
            ys = [float(sp.N(point[1])) for point in polygon]
            ax.fill(xs, ys, alpha=alpha, **plot_kwargs)
        for start, end in data.segments:
            ax.plot(
                [_numeric_pair(start)[0], _numeric_pair(end)[0]],
                [_numeric_pair(start)[1], _numeric_pair(end)[1]],
                **plot_kwargs,
            )
        for point in data.points:
            ax.plot(
                [float(sp.N(point[0]))],
                [float(sp.N(point[1]))],
                marker="o",
                linestyle="None",
                **plot_kwargs,
            )
        ax.set_xlabel(sp.sstr(data.variables[0]))
        ax.set_ylabel(sp.sstr(data.variables[1]))
        ax.set_aspect("equal", adjustable="box")
    else:
        raise NotImplementedError("plot_solution currently supports only 1D and 2D solution views")
    if bounds and len(bounds) >= 2:
        ax.set_xlim(float(sp.N(bounds[0][0])), float(sp.N(bounds[0][1])))
        ax.set_ylim(float(sp.N(bounds[1][0])), float(sp.N(bounds[1][1])))
    if show:
        plt.show()
    return ax


__all__ = ["SolutionPlotData", "discretize_solution", "plot_solution"]


def discretize_region_geometry(
    region, *, variables=None, samples_per_curve: int = 64
) -> SolutionPlotData:
    """Return lightweight geometry for explicit standard-region objects.

    This is a plotting/discretization front end for the explicit region classes
    in :mod:`semialg.standard_regions`. It is intentionally lightweight: it
    returns points, segments, and polygons rather than a certified mesh.
    """

    from .standard_regions import (
        BallRegion,
        BoxRegion,
        CapsuleRegion,
        IntervalRegion,
        ParallelogramRegion,
        PointRegion,
        PolygonRegion,
        SphereRegion,
        StadiumRegion,
        StandardRegion,
    )

    if not isinstance(region, StandardRegion):
        raise TypeError("discretize_region_geometry expects a StandardRegion object")
    n = region.ambient_dimension()
    if variables is None:
        variables = tuple(sp.Symbol(f"x{i + 1}", real=True) for i in range(n))
    else:
        variables = tuple(variables)

    points = []
    segments = []
    polygons = []
    if isinstance(region, PointRegion):
        points.extend(region.points)
    elif isinstance(region, IntervalRegion):
        if region.lower == region.upper:
            points.append((region.lower,))
        else:
            segments.append(((region.lower,), (region.upper,)))
    elif isinstance(region, BoxRegion) and len(region.bounds) == 2:
        (x0, x1), (y0, y1) = region.bounds
        polygons.append(((x0, y0), (x1, y0), (x1, y1), (x0, y1)))
    elif isinstance(region, PolygonRegion):
        polygons.append(region.vertices)
    elif isinstance(region, ParallelogramRegion):
        o = sp.Matrix(region.origin)
        v1 = sp.Matrix(region.vectors[0])
        v2 = sp.Matrix(region.vectors[1])
        verts = [tuple(o), tuple(o + v1), tuple(o + v1 + v2), tuple(o + v2)]
        if len(verts[0]) == 2:
            polygons.append(tuple(verts))
    elif isinstance(region, (BallRegion, SphereRegion)) and len(region.center) == 2:
        cx, cy = region.center
        r = region.radius
        pts = tuple(
            (
                sp.simplify(cx + r * sp.cos(2 * sp.pi * i / samples_per_curve)),
                sp.simplify(cy + r * sp.sin(2 * sp.pi * i / samples_per_curve)),
            )
            for i in range(samples_per_curve)
        )
        if isinstance(region, SphereRegion):
            segments.extend((pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts)))
        else:
            polygons.append(pts)
    elif isinstance(region, (StadiumRegion, CapsuleRegion)) and len(region.start) == 2:
        # Coarse convex-hull style representation; enough for visual debugging.
        x0, y0 = region.start
        x1, y1 = region.end
        r = region.radius
        polygons.append(((x0, y0 - r), (x1, y1 - r), (x1, y1 + r), (x0, y0 + r)))
    else:
        raise NotImplementedError(f"discretization is not implemented for {type(region).__name__}")

    return SolutionPlotData(
        variables=tuple(variables),
        points=tuple(points),
        segments=tuple(segments),
        polygons=tuple(polygons),
        samples=(),
        source=f"standard_region:{type(region).__name__}",
        diagnostics={"samples_per_curve": samples_per_curve},
    )


def plot_region_geometry(
    region,
    *,
    variables=None,
    ax=None,
    show: bool = False,
    samples_per_curve: int = 64,
    **plot_kwargs,
):
    """Plot an explicit standard-region object using Matplotlib."""

    dummy = type("_SolutionLike", (), {})()
    data = discretize_region_geometry(
        region, variables=variables, samples_per_curve=samples_per_curve
    )
    dummy.variables = data.variables
    dummy.points = data.points
    dummy.segments = data.segments
    dummy.polygons = data.polygons
    # Reuse simple Matplotlib logic by drawing directly.
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise ImportError("plot_region_geometry requires matplotlib") from exc
    if ax is None:
        _, ax = plt.subplots()
    if len(data.variables) == 1:
        for point in data.points:
            ax.plot([float(sp.N(point[0]))], [0.0], marker="o", linestyle="None", **plot_kwargs)
        for start, end in data.segments:
            ax.plot([float(sp.N(start[0])), float(sp.N(end[0]))], [0.0, 0.0], **plot_kwargs)
    elif len(data.variables) == 2:
        for polygon in data.polygons:
            xs = [float(sp.N(point[0])) for point in polygon]
            ys = [float(sp.N(point[1])) for point in polygon]
            ax.fill(xs, ys, alpha=plot_kwargs.pop("alpha", 0.25), **plot_kwargs)
        for start, end in data.segments:
            ax.plot(
                [float(sp.N(start[0])), float(sp.N(end[0]))],
                [float(sp.N(start[1])), float(sp.N(end[1]))],
                **plot_kwargs,
            )
        for point in data.points:
            ax.plot(
                [float(sp.N(point[0]))],
                [float(sp.N(point[1]))],
                marker="o",
                linestyle="None",
                **plot_kwargs,
            )
        ax.set_aspect("equal", adjustable="box")
    else:
        raise NotImplementedError("plot_region_geometry currently supports 1D and 2D views")
    if show:
        plt.show()
    return ax


__all__ = [
    "SolutionPlotData",
    "discretize_solution",
    "plot_solution",
    "discretize_region_geometry",
    "plot_region_geometry",
]
