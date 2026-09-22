"""Coercion of explicit/formula regions to :class:`SemialgebraicRegion`.

Kept separate from the symbolic region value so importing the core abstraction
does not eagerly import the full explicit-geometry hierarchy.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from ._zero_testing import certified_zero
from .normalization import normalize_formula, normalize_variables
from .quantifiers import Exists
from .symbolic_regions import SemialgebraicRegion, _default_variables, _fresh_symbols


def _composite_region_formula(region, variables: tuple[sp.Symbol, ...]) -> sp.Expr | None:
    """Lower parametric, transformed, and Boolean composition representations."""
    from .standard_regions import BooleanRegion, ParametricRegion, TransformedRegion

    if isinstance(region, ParametricRegion):
        body = sp.And(
            normalize_formula(region.assumptions),
            *(sp.And(param >= lo, param <= hi) for param, lo, hi in region.limits),
            *(sp.Eq(v, mapping) for v, mapping in zip(variables, region.mapping, strict=True)),
        )
        return Exists(region.parameters, body)
    if isinstance(region, TransformedRegion):
        base_vars = tuple(region.base_variables)
        if len(base_vars) != region.base.ambient_dimension():
            raise ValueError("TransformedRegion base variable count does not match base dimension")
        base_formula = _standard_region_formula(region.base, base_vars)
        body = sp.And(
            base_formula,
            *(sp.Eq(v, mapping) for v, mapping in zip(variables, region.mapping, strict=True)),
        )
        return Exists(base_vars, body)
    if isinstance(region, BooleanRegion):
        formulas = [as_semialgebraic_region(r, variables).formula for r in region.regions]
        if region.op == "union":
            return sp.Or(*formulas)
        if region.op == "intersection":
            return sp.And(*formulas)
        if region.op == "difference":
            return sp.And(formulas[0], sp.Not(formulas[1]))
        if region.op == "symmetric_difference":
            return sp.Xor(formulas[0], formulas[1])
        if region.op == "complement":
            return sp.Not(formulas[0])
    return None


def _swept_region_formula(region, variables: tuple[sp.Symbol, ...]) -> sp.Expr | None:
    """Lower segment-swept and axial solid families without affecting shape dispatch."""
    from .standard_regions import Capsule, Cone, Cylinder, Stadium

    if isinstance(region, (Stadium, Capsule)):
        t = sp.Dummy("t", real=True)
        axis_point = tuple(a + t * (b - a) for a, b in zip(region.start, region.end, strict=True))
        radial = sum((v - q) ** 2 for v, q in zip(variables, axis_point, strict=True))
        return Exists(t, sp.And(t >= 0, t <= 1, radial <= region.radius**2))
    if isinstance(region, (Cylinder, Cone)):
        t = sp.Dummy("t", real=True)
        direction = tuple(b - a for a, b in zip(region.start, region.end, strict=True))
        norm_sq = sp.expand(sum(delta**2 for delta in direction))
        if certified_zero(norm_sq) is True:
            radial = sum((v - a) ** 2 for v, a in zip(variables, region.start, strict=True))
            return radial <= region.radius**2
        offset = tuple(v - a for v, a in zip(variables, region.start, strict=True))
        projection = sp.expand(sum(o * d for o, d in zip(offset, direction, strict=True)))
        perpendicular_sq = sp.expand(sum(o**2 for o in offset) - t**2 * norm_sq)
        radius_sq = (
            region.radius**2
            if isinstance(region, Cylinder) and not isinstance(region, Cone)
            else ((1 - t) * region.radius) ** 2
        )
        return Exists(
            t, sp.And(t >= 0, t <= 1, sp.Eq(projection, t * norm_sq), perpendicular_sq <= radius_sq)
        )
    return None


def _standard_region_formula(region, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    """Lower a ``StandardRegion`` to a symbolic semialgebraic membership formula."""

    from .standard_regions import (
        AffineHalfSpace,
        AffineSpace,
        Ball,
        Box,
        Ellipsoid,
        EllipsoidBoundary,
        FilledTorus,
        FinitePointSet,
        HalfSpace,
        Hyperplane,
        Interval,
        Parallelepiped,
        Parallelogram,
        Polygon,
        PolyhedralCone,
        Polytope,
        Ray,
        Simplex,
        Sphere,
        SphericalShell,
        TetrahedralComplex,
        Torus,
    )

    if len(variables) != region.ambient_dimension():
        raise ValueError("variable count does not match explicit region ambient dimension")

    if isinstance(region, FinitePointSet):
        pieces = [
            sp.And(*(sp.Eq(v, c) for v, c in zip(variables, p, strict=True))) for p in region.points
        ]
        return sp.Or(*pieces) if pieces else sp.false
    if isinstance(region, Hyperplane):
        return sp.Eq(
            sum(
                n * (v - p) for n, v, p in zip(region.normal, variables, region.point, strict=True)
            ),
            0,
        )
    if isinstance(region, HalfSpace):
        return (
            sum(n * (v - p) for n, v, p in zip(region.normal, variables, region.point, strict=True))
            <= 0
        )
    if isinstance(region, Ray):
        t = sp.Dummy("t", real=True)
        return Exists(
            t,
            sp.And(
                t >= 0,
                *(
                    sp.Eq(v, p + t * d)
                    for v, p, d in zip(variables, region.point, region.direction, strict=True)
                ),
            ),
        )
    if isinstance(region, AffineHalfSpace):
        params = _fresh_symbols("u", len(region.directions))
        t = sp.Dummy("t", real=True)
        coords = [
            sp.Eq(
                v,
                region.point[j]
                + sum(a * d[j] for a, d in zip(params, region.directions, strict=True))
                + t * region.inward[j],
            )
            for j, v in enumerate(variables)
        ]
        return Exists((*params, t), sp.And(t >= 0, *coords))
    if isinstance(region, AffineSpace):
        params = _fresh_symbols("u", len(region.directions))
        if not params:
            return sp.And(*(sp.Eq(v, p) for v, p in zip(variables, region.point, strict=True)))
        coords = [
            sp.Eq(
                v,
                region.point[j]
                + sum(a * d[j] for a, d in zip(params, region.directions, strict=True)),
            )
            for j, v in enumerate(variables)
        ]
        return Exists(params, sp.And(*coords))
    if isinstance(region, PolyhedralCone):
        free = _fresh_symbols("u", len(region.directions))
        nonnegative = _fresh_symbols("t", len(region.rays))
        coords = [
            sp.Eq(
                v,
                region.point[j]
                + sum(a * d[j] for a, d in zip(free, region.directions, strict=True))
                + sum(t * ray[j] for t, ray in zip(nonnegative, region.rays, strict=True)),
            )
            for j, v in enumerate(variables)
        ]
        body = sp.And(*(t >= 0 for t in nonnegative), *coords)
        params = (*free, *nonnegative)
        return Exists(params, body) if params else body
    if isinstance(region, Polytope):
        weights = _fresh_symbols("lambda", len(region.vertices))
        coords = [
            sp.Eq(
                v,
                sum(w * vertex[j] for w, vertex in zip(weights, region.vertices, strict=True)),
            )
            for j, v in enumerate(variables)
        ]
        return Exists(
            weights,
            sp.And(sp.Eq(sum(weights), 1), *(w >= 0 for w in weights), *coords),
        )
    if isinstance(region, Interval):
        (x,) = variables
        lo = x >= region.lower if region.lower_closed else x > region.lower
        hi = x <= region.upper if region.upper_closed else x < region.upper
        return sp.And(lo, hi)
    if isinstance(region, Box):
        return sp.And(
            *(
                sp.And(v >= lo, v <= hi)
                for v, (lo, hi) in zip(variables, region.bounds, strict=True)
            )
        )
    if isinstance(region, Sphere):
        radial = sum((v - c) ** 2 for v, c in zip(variables, region.center, strict=True))
        return sp.Eq(radial, region.radius**2)
    if isinstance(region, Ball):
        radial = sum((v - c) ** 2 for v, c in zip(variables, region.center, strict=True))
        return radial <= region.radius**2
    if isinstance(region, (Ellipsoid, EllipsoidBoundary)):
        delta = sp.Matrix(variables) - sp.Matrix(region.center)
        quadratic = sp.expand((delta.T * region.shape_matrix.inv() * delta)[0])
        return sp.Eq(quadratic, 1) if isinstance(region, EllipsoidBoundary) else quadratic <= 1
    if isinstance(region, SphericalShell):
        radial = sum((v - c) ** 2 for v, c in zip(variables, region.center, strict=True))
        return sp.And(radial >= region.inner_radius**2, radial <= region.outer_radius**2)
    if isinstance(region, Torus):
        x, y, z = variables
        c1, c2, c3 = region.center
        rho_sq = (x - c1) ** 2 + (y - c2) ** 2
        z_sq = (z - c3) ** 2
        r_major = region.major_radius
        r_minor = region.minor_radius
        quartic = sp.expand(
            (rho_sq + z_sq + r_major**2 - r_minor**2) ** 2 - 4 * r_major**2 * rho_sq
        )
        return quartic <= 0 if isinstance(region, FilledTorus) else sp.Eq(quartic, 0)
    if isinstance(region, Simplex):
        lambdas = _fresh_symbols("lambda", len(region.vertices))
        coords = [
            sp.Eq(
                v,
                sum(lam * vertex[j] for lam, vertex in zip(lambdas, region.vertices, strict=True)),
            )
            for j, v in enumerate(variables)
        ]
        body = sp.And(
            sp.Eq(sum(lambdas), 1),
            *(lam >= 0 for lam in lambdas),
            *coords,
        )
        return Exists(lambdas, body)
    if isinstance(region, Polygon):
        return sp.Or(
            *(_standard_region_formula(triangle, variables) for triangle in region.triangulation())
        )
    if isinstance(region, TetrahedralComplex):
        return sp.Or(*(_standard_region_formula(tet, variables) for tet in region.tetrahedra))
    if isinstance(region, (Parallelogram, Parallelepiped)):
        params = _fresh_symbols("u", len(region.vectors))
        body = sp.And(
            *(sp.And(param >= 0, param <= 1) for param in params),
            *(
                sp.Eq(
                    v,
                    region.origin[j]
                    + sum(
                        param * vec[j] for param, vec in zip(params, region.vectors, strict=True)
                    ),
                )
                for j, v in enumerate(variables)
            ),
        )
        return Exists(params, body)
    composite = _composite_region_formula(region, variables)
    if composite is not None:
        return composite
    swept = _swept_region_formula(region, variables)
    if swept is not None:
        return swept
    raise NotImplementedError(f"cannot lower {type(region).__name__} to a semialgebraic formula")


def as_semialgebraic_region(
    region: object,
    variables: Sequence[sp.Symbol | str] | None = None,
) -> SemialgebraicRegion:
    """Coerce a formula or explicit :class:`Geometry` to ``SemialgebraicRegion``."""

    if isinstance(region, SemialgebraicRegion):
        if variables is None:
            return region
        vars_ = tuple(normalize_variables(variables, region.formula, append_context_symbols=False))
        if vars_ == region.variables:
            return region
        if len(vars_) != len(region.variables):
            raise ValueError("replacement variable count does not match region dimension")
        formula = region.formula.xreplace(dict(zip(region.variables, vars_, strict=True)))
        return SemialgebraicRegion(formula, vars_)

    from .standard_regions import Geometry

    if isinstance(region, Geometry):
        vars_ = (
            tuple(normalize_variables(variables, append_context_symbols=False))
            if variables is not None
            else _default_variables(region.ambient_dimension())
        )
        if len(vars_) != region.ambient_dimension():
            raise ValueError("variable count does not match explicit region ambient dimension")
        return SemialgebraicRegion(_standard_region_formula(region, vars_), vars_)

    expr = normalize_formula(region)
    vars_ = tuple(normalize_variables(variables, expr, append_context_symbols=variables is None))
    return SemialgebraicRegion(expr, vars_)


__all__ = ["as_semialgebraic_region"]


# ``symbolic_regions`` historically re-exports this operation.  Install the
# canonical implementation object after this module is fully initialized so the
# re-export preserves object identity without creating a second implementation
# or an eager circular import.
from . import symbolic_regions as _symbolic_regions  # noqa: E402

_symbolic_regions.as_semialgebraic_region = as_semialgebraic_region
