"""Exact convex hulls of finite point sets."""

from __future__ import annotations

from functools import cmp_to_key

import sympy as sp

from ._zero_testing import certified_zero


def _unique_points(points):
    pts = tuple(tuple(sp.sympify(x) for x in p) for p in points)
    if not pts:
        raise ValueError("convex_hull requires at least one point")
    n = len(pts[0])
    if any(len(p) != n for p in pts):
        raise ValueError("all points must have the same ambient dimension")
    return tuple(dict.fromkeys(pts))


def _intrinsic_chart(points):
    anchor = sp.Matrix(points[0])
    if len(points) == 1:
        return anchor, sp.zeros(len(anchor), 0), ((),)
    D = sp.Matrix.hstack(*(sp.Matrix(p) - anchor for p in points[1:]))
    basis_cols = D.columnspace()
    B = sp.Matrix.hstack(*basis_cols) if basis_cols else sp.zeros(len(anchor), 0)
    d = B.cols
    if d == 0:
        return anchor, B, tuple(() for _ in points)
    pivot_rows = tuple(sp.Matrix(B.T).rref()[1])
    square = B[list(pivot_rows), :]
    if certified_zero(square.det()) is not False:
        raise NotImplementedError("could not certify an intrinsic affine chart")
    inverse = square.inv()
    coords = []
    for p in points:
        delta = sp.Matrix(p) - anchor
        q = inverse * delta[list(pivot_rows), :]
        if any(certified_zero(v) is not True for v in (B * q - delta)):
            raise ValueError("points do not lie in the computed affine hull")
        coords.append(tuple(sp.simplify(q[i, 0]) for i in range(d)))
    return anchor, B, tuple(coords)


def convex_hull(points, *, canonical=True):
    """Return the exact convex hull of a finite point set.

    Redundant and interior input points are removed exactly.  Lower-dimensional
    point sets are solved in an exact intrinsic affine chart and mapped back.
    """
    from .polyhedral import h_representation_from_vertices
    from .standard_regions import FinitePointSet, Polytope

    pts = _unique_points(points)
    anchor, basis, coords = _intrinsic_chart(pts)
    d = basis.cols
    if d == 0:
        return FinitePointSet((pts[0],))
    if d == 1:
        from .exact_arithmetic import compare_exact_reals

        ordered = sorted(
            coords,
            key=cmp_to_key(lambda a, b: compare_exact_reals(a[0], b[0])),
        )
        extreme_coords = (ordered[0], ordered[-1])
    else:
        hrep = h_representation_from_vertices(coords)
        extreme_coords = hrep.vertices()
    vertices = []
    for q in extreme_coords:
        p = anchor + basis * sp.Matrix(q)
        point = tuple(sp.simplify(p[i, 0]) for i in range(p.rows))
        if point not in vertices:
            vertices.append(point)
    polytope = Polytope(vertices)
    if not canonical:
        return polytope
    from .canonicalization import canonicalize_region

    return canonicalize_region(polytope)


__all__ = ["convex_hull"]
