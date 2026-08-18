from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from .parametric_integration import integrate_over_parametric_region
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
    StandardRegion,
    TetrahedronRegion,
    TransformedRegion,
)


def _symbols(names: Sequence[sp.Symbol | str]) -> tuple[sp.Symbol, ...]:
    return tuple(sp.Symbol(v, real=True) if isinstance(v, str) else v for v in names)


def _monomial_integral_box(
    exponents: tuple[int, ...], bounds: tuple[tuple[sp.Expr, sp.Expr], ...]
) -> sp.Expr:
    value = sp.Integer(1)
    for exp, (lo, hi) in zip(exponents, bounds, strict=True):
        value *= sp.simplify((hi ** (exp + 1) - lo ** (exp + 1)) / (exp + 1))
    return sp.simplify(value)


def _integrate_polynomial_over_box(
    expr: sp.Expr, variables: tuple[sp.Symbol, ...], bounds: tuple[tuple[sp.Expr, sp.Expr], ...]
) -> sp.Expr | None:
    try:
        poly = sp.Poly(sp.expand(expr), *variables)
    except Exception:
        return None
    total = sp.Integer(0)
    for monom, coeff in poly.terms():
        total += coeff * _monomial_integral_box(tuple(monom), bounds)
    return sp.simplify(total)


def _simplex_parametric_region(region: SimplexRegion) -> ParametricRegion:
    vertices = region.vertices
    k = len(vertices) - 1
    if k == 0:
        return ParametricRegion((), (), vertices[0])
    params = tuple(sp.Symbol(f"_u{i + 1}", real=True) for i in range(k))
    p0 = sp.Matrix(vertices[0])
    mapping_vec = p0
    for i, u in enumerate(params):
        mapping_vec += u * (sp.Matrix(vertices[i + 1]) - p0)
    limits = []
    # u_k inner first: 0 <= u_k <= 1-u_1-...-u_{k-1}
    for i in reversed(range(k)):
        upper = 1 - sum(params[:i])
        limits.append((params[i], sp.Integer(0), upper))
    return ParametricRegion(params, tuple(limits), tuple(mapping_vec))


def _parallelepiped_parametric_region(
    region: ParallelepipedRegion | ParallelogramRegion,
) -> ParametricRegion:
    vectors = tuple(region.vectors)  # type: ignore[attr-defined]
    params = tuple(sp.Symbol(f"_u{i + 1}", real=True) for i in range(len(vectors)))
    origin = sp.Matrix(region.origin)  # type: ignore[attr-defined]
    mapping = origin
    for u, vec in zip(params, vectors, strict=True):
        mapping += u * sp.Matrix(vec)
    limits = tuple((u, sp.Integer(0), sp.Integer(1)) for u in reversed(params))
    return ParametricRegion(params, limits, tuple(mapping))


def _ball_monomial_centered_integral(exponents: tuple[int, ...], radius: sp.Expr) -> sp.Expr:
    if any(e % 2 for e in exponents):
        return sp.Integer(0)
    n = len(exponents)
    degree = sum(exponents)
    angular = sp.prod(sp.gamma(sp.Rational(e + 1, 2)) for e in exponents)
    return sp.simplify(radius ** (degree + n) * angular / sp.gamma(1 + sp.Rational(degree + n, 2)))


def _sphere_monomial_centered_integral(exponents: tuple[int, ...], radius: sp.Expr) -> sp.Expr:
    if any(e % 2 for e in exponents):
        return sp.Integer(0)
    n = len(exponents)
    degree = sum(exponents)
    angular = 2 * sp.prod(sp.gamma(sp.Rational(e + 1, 2)) for e in exponents)
    return sp.simplify(radius ** (degree + n - 1) * angular / sp.gamma(sp.Rational(degree + n, 2)))


def _integrate_polynomial_over_centered_region(
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    center: tuple[sp.Expr, ...],
    radius: sp.Expr,
    *,
    surface: bool = False,
) -> sp.Expr | None:
    shifted_vars = tuple(sp.Symbol(f"_{v.name}_centered", real=True) for v in variables)
    shifted = sp.expand(
        expr.subs({v: c + u for v, c, u in zip(variables, center, shifted_vars, strict=True)})
    )
    try:
        poly = sp.Poly(shifted, *shifted_vars)
    except Exception:
        return None
    total = sp.Integer(0)
    for monom, coeff in poly.terms():
        if surface:
            total += coeff * _sphere_monomial_centered_integral(tuple(monom), radius)
        else:
            total += coeff * _ball_monomial_centered_integral(tuple(monom), radius)
    return sp.simplify(total)


def _segment_length(start: tuple[sp.Expr, ...], end: tuple[sp.Expr, ...]) -> sp.Expr:
    return sp.sqrt(sum((b - a) ** 2 for a, b in zip(start, end, strict=True)))


def _aligned_z_parametric_cylinder(
    region: CylinderRegion, cone: bool = False
) -> ParametricRegion | None:
    if len(region.start) != 3 or len(region.end) != 3:
        return None
    x0, y0, z0 = region.start
    x1, y1, z1 = region.end
    if sp.simplify(x1 - x0) != 0 or sp.simplify(y1 - y0) != 0:
        return None
    r, th, h = sp.symbols("_r _theta _h", real=True)
    height = sp.simplify(z1 - z0)
    if cone:
        rr = region.radius * (1 - h)
    else:
        rr = region.radius
    mapping = (x0 + rr * r * sp.cos(th), y0 + rr * r * sp.sin(th), z0 + height * h)
    # Mapping is not algebraic, but parametric metric integration can handle it for many polynomials.
    return ParametricRegion((r, th, h), ((r, 0, 1), (th, 0, 2 * sp.pi), (h, 0, 1)), mapping)


def _integrate_boolean_region(
    integrand: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    region: BooleanRegion,
    *,
    method: str,
    precision: int,
) -> sp.Expr:
    if region.op == "union":
        if region.assume_disjoint:
            return sp.simplify(
                sum(
                    integrate_over_standard_region(
                        integrand, r, variables, method=method, precision=precision
                    )
                    for r in region.regions
                )
            )
        if len(region.regions) == 2:
            a, b = region.regions
            return sp.simplify(
                integrate_over_standard_region(
                    integrand, a, variables, method=method, precision=precision
                )
                + integrate_over_standard_region(
                    integrand, b, variables, method=method, precision=precision
                )
                - integrate_over_standard_region(
                    integrand,
                    BooleanRegion("intersection", (a, b)),
                    variables,
                    method=method,
                    precision=precision,
                )
            )
        raise NotImplementedError(
            "general non-disjoint unions require explicit disjoint decomposition"
        )
    if region.op == "difference":
        if len(region.regions) != 2:
            raise ValueError("difference requires two regions")
        a, b = region.regions
        return sp.simplify(
            integrate_over_standard_region(
                integrand, a, variables, method=method, precision=precision
            )
            - integrate_over_standard_region(
                integrand,
                BooleanRegion("intersection", (a, b)),
                variables,
                method=method,
                precision=precision,
            )
        )
    if region.op == "symmetric_difference":
        if len(region.regions) != 2:
            raise ValueError("symmetric_difference requires two regions")
        a, b = region.regions
        return sp.simplify(
            integrate_over_standard_region(
                integrand, a, variables, method=method, precision=precision
            )
            + integrate_over_standard_region(
                integrand, b, variables, method=method, precision=precision
            )
            - 2
            * integrate_over_standard_region(
                integrand,
                BooleanRegion("intersection", (a, b)),
                variables,
                method=method,
                precision=precision,
            )
        )
    if region.op == "intersection":
        # First useful exact intersection: intervals and boxes.
        if all(isinstance(r, IntervalRegion) for r in region.regions) and len(variables) == 1:
            lo = max((r.lower for r in region.regions), key=lambda z: float(sp.N(z)))
            hi = min((r.upper for r in region.regions), key=lambda z: float(sp.N(z)))
            if bool(sp.simplify(hi - lo) < 0):
                return sp.Integer(0)
            return integrate_over_standard_region(
                integrand, IntervalRegion(lo, hi), variables, method=method, precision=precision
            )
        if all(isinstance(r, BoxRegion) for r in region.regions):
            bounds = []
            for i in range(len(variables)):
                lo = max((r.bounds[i][0] for r in region.regions), key=lambda z: float(sp.N(z)))
                hi = min((r.bounds[i][1] for r in region.regions), key=lambda z: float(sp.N(z)))
                if bool(sp.simplify(hi - lo) < 0):
                    return sp.Integer(0)
                bounds.append((lo, hi))
            return integrate_over_standard_region(
                integrand, BoxRegion(bounds), variables, method=method, precision=precision
            )
        raise NotImplementedError("this BooleanRegion intersection is not yet supported exactly")
    raise NotImplementedError(f"unsupported BooleanRegion op: {region.op}")


def integrate_over_standard_region(
    integrand: object,
    region: StandardRegion,
    variables: Sequence[sp.Symbol | str],
    *,
    method: str = "symbolic",
    precision: int = 50,
) -> sp.Expr:
    """Integrate over an explicit standard region object.

    This layer complements semialgebraic-formula integration with exact formulas
    and parametrizations for common geometric regions.
    """

    vars_ = _symbols(variables)
    expr = sp.sympify(integrand)

    if isinstance(region, PointRegion):
        return sp.simplify(sum(expr.subs(dict(zip(vars_, p, strict=True))) for p in region.points))
    if isinstance(region, IntervalRegion):
        value = sp.integrate(expr, (vars_[0], region.lower, region.upper))
        if isinstance(value, sp.Integral) or value.has(sp.Integral):
            if method == "symbolic":
                raise NotImplementedError("could not integrate over interval symbolically")
            return value.evalf(precision)
        return sp.simplify(value)
    if isinstance(region, BoxRegion):
        fast = _integrate_polynomial_over_box(expr, vars_, region.bounds)
        if fast is not None:
            return fast
        value = sp.integrate(
            expr,
            *tuple(
                (v, lo, hi)
                for v, (lo, hi) in reversed(tuple(zip(vars_, region.bounds, strict=True)))
            ),
        )
        return sp.N(value, precision) if method == "numeric" else sp.simplify(value)
    if isinstance(region, (SimplexRegion, TetrahedronRegion)):
        pregion = _simplex_parametric_region(region)
        return integrate_over_parametric_region(
            expr, vars_, pregion, method=method, precision=precision
        )  # type: ignore[return-value]
    if isinstance(region, PolygonRegion):
        return sp.simplify(
            sum(
                integrate_over_standard_region(expr, tri, vars_, method=method, precision=precision)
                for tri in region.triangulation()
            )
        )
    if isinstance(region, PolyhedronRegion):
        return sp.simplify(
            sum(
                integrate_over_standard_region(expr, tet, vars_, method=method, precision=precision)
                for tet in region.tetrahedra
            )
        )
    if isinstance(region, (ParallelogramRegion, ParallelepipedRegion)):
        return integrate_over_parametric_region(
            expr,
            vars_,
            _parallelepiped_parametric_region(region),
            method=method,
            precision=precision,
        )  # type: ignore[return-value]
    if isinstance(region, PrismRegion):
        base_pieces = (
            region.base.triangulation()
            if isinstance(region.base, PolygonRegion)
            else (region.base,)
        )
        total = sp.Integer(0)
        for tri in base_pieces:
            # triangular prism is parametrized as triangle + w*vector
            simplex_param = _simplex_parametric_region(tri)
            w = sp.Symbol("_w", real=True)
            mapping = tuple(
                sp.simplify(c + w * dv)
                for c, dv in zip(simplex_param.mapping, region.vector, strict=True)
            )
            pregion = ParametricRegion(
                (*simplex_param.parameters, w), (*simplex_param.limits, (w, 0, 1)), mapping
            )
            total += integrate_over_parametric_region(
                expr, vars_, pregion, method=method, precision=precision
            )  # type: ignore[operator]
        return sp.simplify(total)
    if isinstance(region, PyramidRegion):
        base_pieces = (
            region.base.triangulation()
            if isinstance(region.base, PolygonRegion)
            else (region.base,)
        )
        total = sp.Integer(0)
        for tri in base_pieces:
            tet = TetrahedronRegion((*tri.vertices, region.apex))
            total += integrate_over_standard_region(
                expr, tet, vars_, method=method, precision=precision
            )
        return sp.simplify(total)
    if isinstance(region, BallRegion) and not isinstance(region, SphereRegion):
        fast = _integrate_polynomial_over_centered_region(
            expr, vars_, region.center, region.radius, surface=False
        )
        if fast is not None:
            return fast
    if isinstance(region, SphereRegion):
        fast = _integrate_polynomial_over_centered_region(
            expr, vars_, region.center, region.radius, surface=True
        )
        if fast is not None:
            return fast
    if isinstance(region, SphericalShellRegion):
        outer = _integrate_polynomial_over_centered_region(
            expr, vars_, region.center, region.outer_radius, surface=False
        )
        inner = _integrate_polynomial_over_centered_region(
            expr, vars_, region.center, region.inner_radius, surface=False
        )
        if outer is not None and inner is not None:
            return sp.simplify(outer - inner)
    if isinstance(region, CylinderRegion) and not isinstance(region, ConeRegion):
        if expr == 1:
            return sp.simplify(sp.pi * region.radius**2 * _segment_length(region.start, region.end))
        pregion = _aligned_z_parametric_cylinder(region, cone=False)
        if pregion is not None:
            return integrate_over_parametric_region(
                expr, vars_, pregion, method=method, precision=precision
            )  # type: ignore[return-value]
    if isinstance(region, ConeRegion):
        if expr == 1:
            return sp.simplify(
                sp.pi * region.radius**2 * _segment_length(region.start, region.end) / 3
            )
        pregion = _aligned_z_parametric_cylinder(region, cone=True)
        if pregion is not None:
            return integrate_over_parametric_region(
                expr, vars_, pregion, method=method, precision=precision
            )  # type: ignore[return-value]
    if isinstance(region, StadiumRegion) and not isinstance(region, CapsuleRegion):
        if expr == 1:
            return sp.simplify(
                2 * region.radius * _segment_length(region.start, region.end)
                + sp.pi * region.radius**2
            )
    if isinstance(region, CapsuleRegion):
        if expr == 1:
            n = len(region.start)
            length = _segment_length(region.start, region.end)
            if n == 2:
                return sp.simplify(2 * region.radius * length + sp.pi * region.radius**2)
            if n == 3:
                return sp.simplify(
                    sp.pi * region.radius**2 * length + sp.Rational(4, 3) * sp.pi * region.radius**3
                )
    if isinstance(region, ParametricRegion):
        return integrate_over_parametric_region(
            expr, vars_, region, method=method, precision=precision
        )  # type: ignore[return-value]
    if isinstance(region, TransformedRegion):
        # Transform by composing an available parametrization of the base when possible.
        base_vars = region.base_variables
        if isinstance(region.base, BoxRegion):
            params = base_vars
            limits = tuple(
                (v, lo, hi)
                for v, (lo, hi) in reversed(tuple(zip(params, region.base.bounds, strict=True)))
            )
            pregion = ParametricRegion(params, limits, region.mapping)
            return integrate_over_parametric_region(
                expr, vars_, pregion, method=method, precision=precision
            )  # type: ignore[return-value]
        if isinstance(region.base, (SimplexRegion, TetrahedronRegion)):
            base_param = _simplex_parametric_region(region.base)
            composed = tuple(
                sp.simplify(e.subs(dict(zip(base_vars, base_param.mapping, strict=True))))
                for e in region.mapping
            )
            return integrate_over_parametric_region(
                expr,
                vars_,
                ParametricRegion(base_param.parameters, base_param.limits, composed),
                method=method,
                precision=precision,
            )  # type: ignore[return-value]
        raise NotImplementedError(
            "TransformedRegion currently supports BoxRegion and SimplexRegion bases"
        )
    if isinstance(region, BooleanRegion):
        return _integrate_boolean_region(expr, vars_, region, method=method, precision=precision)

    raise NotImplementedError(
        f"standard region integration is not implemented for {type(region).__name__}"
    )


__all__ = ["integrate_over_standard_region"]
