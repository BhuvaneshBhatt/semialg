"""Coercion of explicit/formula regions to :class:`SemialgebraicRegion`.

Kept separate from the symbolic region value so importing the core abstraction
does not eagerly import the full ``StandardRegion`` hierarchy.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from .normalization import normalize_formula, normalize_variables
from .quantifiers import Exists
from .symbolic_regions import SemialgebraicRegion, _default_variables, _fresh_symbols


def _standard_region_formula(region, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    """Lower a ``StandardRegion`` to a symbolic semialgebraic membership formula."""

    from .standard_regions import (
        BallRegion,
        BooleanRegion,
        BoxRegion,
        CapsuleRegion,
        ConeRegion,
        CylinderRegion,
        IntervalRegion,
        ParallelepipedRegion,
        ParallelogramRegion,
        ParametricRegion,
        PointRegion,
        PolygonRegion,
        PolyhedronRegion,
        PrismRegion,
        PyramidRegion,
        SimplexRegion,
        SphereRegion,
        SphericalShellRegion,
        StadiumRegion,
        TransformedRegion,
    )

    if len(variables) != region.ambient_dimension():
        raise ValueError("variable count does not match explicit region ambient dimension")

    if isinstance(region, PointRegion):
        pieces = [
            sp.And(*(sp.Eq(v, c) for v, c in zip(variables, p, strict=True))) for p in region.points
        ]
        return sp.Or(*pieces) if pieces else sp.false
    if isinstance(region, IntervalRegion):
        (x,) = variables
        lo = x >= region.lower if region.lower_closed else x > region.lower
        hi = x <= region.upper if region.upper_closed else x < region.upper
        return sp.And(lo, hi)
    if isinstance(region, BoxRegion):
        return sp.And(
            *(
                sp.And(v >= lo, v <= hi)
                for v, (lo, hi) in zip(variables, region.bounds, strict=True)
            )
        )
    if isinstance(region, SphereRegion):
        radial = sum((v - c) ** 2 for v, c in zip(variables, region.center, strict=True))
        return sp.Eq(radial, region.radius**2)
    if isinstance(region, BallRegion):
        radial = sum((v - c) ** 2 for v, c in zip(variables, region.center, strict=True))
        return radial <= region.radius**2
    if isinstance(region, SphericalShellRegion):
        radial = sum((v - c) ** 2 for v, c in zip(variables, region.center, strict=True))
        return sp.And(radial >= region.inner_radius**2, radial <= region.outer_radius**2)
    if isinstance(region, SimplexRegion):
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
    if isinstance(region, PolygonRegion):
        return sp.Or(
            *(_standard_region_formula(triangle, variables) for triangle in region.triangulation())
        )
    if isinstance(region, PolyhedronRegion):
        return sp.Or(*(_standard_region_formula(tet, variables) for tet in region.tetrahedra))
    if isinstance(region, (ParallelogramRegion, ParallelepipedRegion)):
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
    if isinstance(region, PrismRegion):
        base_vars = _fresh_symbols("b", len(variables))
        t = sp.Dummy("t", real=True)
        base_formula = _standard_region_formula(region.base, tuple(base_vars))
        body = sp.And(
            base_formula,
            t >= 0,
            t <= 1,
            *(
                sp.Eq(v, b + t * delta)
                for v, b, delta in zip(variables, base_vars, region.vector, strict=True)
            ),
        )
        return Exists((*base_vars, t), body)
    if isinstance(region, PyramidRegion):
        base_vars = _fresh_symbols("b", len(variables))
        t = sp.Dummy("t", real=True)
        base_formula = _standard_region_formula(region.base, tuple(base_vars))
        body = sp.And(
            base_formula,
            t >= 0,
            t <= 1,
            *(
                sp.Eq(v, (1 - t) * b + t * apex)
                for v, b, apex in zip(variables, base_vars, region.apex, strict=True)
            ),
        )
        return Exists((*base_vars, t), body)
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
        subregions = [as_semialgebraic_region(r, variables) for r in region.regions]
        formulas = [sub.formula for sub in subregions]
        if region.op == "union":
            return sp.Or(*formulas)
        if region.op == "intersection":
            return sp.And(*formulas)
        if region.op == "difference":
            if len(formulas) != 2:
                raise ValueError("difference requires two regions")
            return sp.And(formulas[0], sp.Not(formulas[1]))
        if region.op == "symmetric_difference":
            if len(formulas) != 2:
                raise ValueError("symmetric difference requires two regions")
            return sp.Xor(formulas[0], formulas[1])
        if region.op == "complement":
            if len(formulas) != 1:
                raise ValueError("complement requires one region")
            return sp.Not(formulas[0])
    if isinstance(region, (StadiumRegion, CapsuleRegion)):
        # Distance-to-segment representation; exact and valid in any ambient dimension.
        t = sp.Dummy("t", real=True)
        axis_point = tuple(a + t * (b - a) for a, b in zip(region.start, region.end, strict=True))
        radial = sum((v - q) ** 2 for v, q in zip(variables, axis_point, strict=True))
        return Exists(t, sp.And(t >= 0, t <= 1, radial <= region.radius**2))
    if isinstance(region, (CylinderRegion, ConeRegion)):
        # Flat-ended right cylinder/cone with the stored start/end as axis.
        t = sp.Dummy("t", real=True)
        direction = tuple(b - a for a, b in zip(region.start, region.end, strict=True))
        norm_sq = sp.expand(sum(delta**2 for delta in direction))
        if sp.simplify(norm_sq) == 0:
            radial = sum((v - a) ** 2 for v, a in zip(variables, region.start, strict=True))
            return radial <= region.radius**2
        offset = tuple(v - a for v, a in zip(variables, region.start, strict=True))
        projection = sp.expand(sum(o * d for o, d in zip(offset, direction, strict=True)))
        perpendicular_sq = sp.expand(sum(o**2 for o in offset) - t**2 * norm_sq)
        radius_sq = (
            region.radius**2
            if isinstance(region, CylinderRegion) and not isinstance(region, ConeRegion)
            else (t * region.radius) ** 2
        )
        return Exists(
            t,
            sp.And(
                t >= 0,
                t <= 1,
                sp.Eq(projection, t * norm_sq),
                perpendicular_sq <= radius_sq,
            ),
        )
    raise NotImplementedError(f"cannot lower {type(region).__name__} to a semialgebraic formula")


def as_semialgebraic_region(
    region: object,
    variables: Sequence[sp.Symbol | str] | None = None,
) -> SemialgebraicRegion:
    """Coerce a formula or explicit ``StandardRegion`` to ``SemialgebraicRegion``."""

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

    from .standard_regions import StandardRegion

    if isinstance(region, StandardRegion):
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
