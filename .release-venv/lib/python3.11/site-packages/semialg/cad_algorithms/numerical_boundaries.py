"""Numerical evaluation of delineable algebraic CAD boundaries."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import sympy as sp


def _numeric_real_roots(
    poly_expr: sp.Expr, var: sp.Symbol, subs: Mapping[sp.Symbol, object], *, precision: int
) -> tuple[float, ...]:
    expr = sp.expand(poly_expr.subs(subs))
    poly = sp.Poly(expr, var)
    roots = poly.nroots(n=precision, maxsteps=200)
    tol = 10.0 ** (-(max(8, precision // 2)))
    real = sorted(float(sp.re(r)) for r in roots if abs(float(sp.im(r))) <= tol)
    return tuple(real)


def evaluate_delineable_curve(
    z_polynomials: Sequence[tuple[sp.Expr, int]],
    y_boundary: tuple[sp.Expr, int],
    variables: Sequence[sp.Symbol],
    x_values: Iterable[object],
    *,
    z_range: tuple[float, float] | None = None,
    precision: int = 30,
) -> tuple[tuple[float, float, float], ...]:
    """Numerically evaluate delineable algebraic space curves.

    ``y_boundary=(g,j)`` denotes the 1-based j-th real root of ``g(x,y)``;
    each ``(f,i)`` denotes the i-th real root in z of ``f(x,y,z)``.
    """
    if len(variables) != 3:
        raise ValueError("curve evaluation requires variables (x, y, z)")
    x, y, z = variables
    g, j = y_boundary
    if j < 1 or any(i < 1 for _, i in z_polynomials):
        raise ValueError("root indices are 1-based positive integers")
    out: list[tuple[float, float, float]] = []
    for xv0 in x_values:
        xv = float(sp.N(xv0, precision))
        yroots = _numeric_real_roots(g, y, {x: xv}, precision=precision)
        if j > len(yroots):
            raise ValueError(
                f"requested y root {j}, but only {len(yroots)} real roots exist at x={xv}"
            )
        yv = yroots[j - 1]
        for f, i in z_polynomials:
            zroots = _numeric_real_roots(f, z, {x: xv, y: yv}, precision=precision)
            if i > len(zroots):
                raise ValueError(f"requested z root {i}, but only {len(zroots)} real roots exist")
            zv = zroots[i - 1]
            if z_range is not None and not (z_range[0] <= zv <= z_range[1]):
                raise ValueError(f"computed z={zv} lies outside requested range {z_range}")
            out.append((xv, yv, zv))
    return tuple(out)


def evaluate_delineable_surfaces(
    z_boundaries: Sequence[tuple[sp.Expr, int]],
    variables: Sequence[sp.Symbol],
    xy_points: Iterable[Sequence[object]],
    *,
    z_range: tuple[float, float] | None = None,
    merge_tolerance: float | None = None,
    precision: int = 30,
) -> tuple[tuple[tuple[float, float, float], ...], tuple[tuple[int, ...], ...]]:
    """Evaluate algebraic z-root surfaces over supplied xy points.

    Returns ``(coordinates, indices)`` where each row of ``indices`` maps one
    boundary polynomial to coordinate indices for all xy points.
    """
    if len(variables) != 3:
        raise ValueError("surface evaluation requires variables (x, y, z)")
    x, y, z = variables
    pts = [(float(sp.N(p[0], precision)), float(sp.N(p[1], precision))) for p in xy_points]
    coords: list[tuple[float, float, float]] = []
    rows: list[list[int]] = [[] for _ in z_boundaries]
    for xv, yv in pts:
        local: list[tuple[int, float]] = []
        for k, (f, i) in enumerate(z_boundaries):
            if i < 1:
                raise ValueError("root indices are 1-based positive integers")
            roots = _numeric_real_roots(f, z, {x: xv, y: yv}, precision=precision)
            if i > len(roots):
                raise ValueError(f"requested z root {i}, but only {len(roots)} real roots exist")
            zv = roots[i - 1]
            if z_range is not None and not (z_range[0] <= zv <= z_range[1]):
                raise ValueError(f"computed z={zv} lies outside requested range {z_range}")
            index = -1
            if merge_tolerance is not None:
                for old_index, old_z in local:
                    if abs(old_z - zv) <= merge_tolerance:
                        index = old_index
                        break
            if index < 0:
                index = len(coords)
                coords.append((xv, yv, zv))
                local.append((index, zv))
            rows[k].append(index)
    return tuple(coords), tuple(tuple(row) for row in rows)


__all__ = ["evaluate_delineable_curve", "evaluate_delineable_surfaces"]
