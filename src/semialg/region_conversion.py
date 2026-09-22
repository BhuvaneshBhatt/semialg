"""Explicit conversions among exact and numerical region representations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from ._zero_testing import certified_equal, certified_sign
from .boundary_topology import PolygonalSet, PolyhedralComponent, PolyhedralShell, Polyhedron
from .canonicalization import canonicalize_region
from .polyhedron_formula import polyhedron_formula
from .standard_regions import Geometry, Polygon, Polytope


@dataclass(frozen=True)
class RegionConversion:
    value: Any
    representation: str
    exact: bool
    certified: bool
    lossy: bool
    source: str


def _orient2(a, b, p):
    return sp.expand((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]))


def _triangle2_formula(vertices, variables):
    a, b, c = vertices
    p = variables
    orient = _orient2(a, b, c)
    tests = [_orient2(a, b, p), _orient2(b, c, p), _orient2(c, a, p)]
    sign = certified_sign(orient)
    if sign == 1:
        return sp.And(*(t >= 0 for t in tests))
    if sign == -1:
        return sp.And(*(t <= 0 for t in tests))
    raise ValueError("triangle orientation is not certified")


def _polygon_formula(poly, variables):
    return sp.Or(*(_triangle2_formula(t.vertices, variables) for t in poly.triangulation()))


def _polygon_boundary_formula(poly, variables):
    x, y = variables
    p = (x, y)
    pieces = []
    vs = poly.vertices
    for a, b in zip(vs, vs[1:] + vs[:1], strict=True):
        col = sp.Eq(_orient2(a, b, p), 0)
        dot = sp.expand((x - a[0]) * (x - b[0]) + (y - a[1]) * (y - b[1]))
        pieces.append(sp.And(col, dot <= 0))
    return sp.Or(*pieces)


def _polygonal_formula(region, variables):
    pieces = []
    for comp in region.components:
        outer = _polygon_formula(comp.outer_boundary, variables)
        holes = []
        for h in comp.hole_boundaries:
            hf = _polygon_formula(h, variables)
            hb = _polygon_boundary_formula(h, variables)
            holes.append(sp.Or(sp.Not(hf), hb))
        pieces.append(sp.And(outer, *holes))
    return sp.Or(*pieces) if pieces else sp.false


def _semialgebraic(region, variables=None):
    from .region_coercion import as_semialgebraic_region

    if isinstance(region, PolygonalSet):
        vars_ = tuple(variables or sp.symbols("x0:2", real=True))
        from .symbolic_regions import SemialgebraicRegion

        return SemialgebraicRegion(_polygonal_formula(region, vars_), vars_)
    if isinstance(region, Polyhedron):
        vars_ = tuple(variables or sp.symbols("x0:3", real=True))
        from .symbolic_regions import SemialgebraicRegion

        return SemialgebraicRegion(polyhedron_formula(region, vars_), vars_)
    return as_semialgebraic_region(region, variables)


def _ordered_facet(ids, verts):
    ids = list(ids)
    if len(ids) <= 3:
        return ids
    # Project a convex planar facet to a coordinate plane on which it remains
    # two-dimensional, then use an exact monotone-chain boundary order.
    a = sp.Matrix(verts[ids[0]])
    normal = None
    for i in range(1, len(ids)):
        for j in range(i + 1, len(ids)):
            n = (sp.Matrix(verts[ids[i]]) - a).cross(sp.Matrix(verts[ids[j]]) - a)
            if any(certified_equal(v, 0) is False for v in n):
                normal = n
                break
        if normal is not None:
            break
    if normal is None:
        raise ValueError("degenerate facet")
    drop = next(k for k in range(3) if certified_equal(normal[k], 0) is False)
    keep = [k for k in range(3) if k != drop]
    pts = {i: (verts[i][keep[0]], verts[i][keep[1]]) for i in ids}
    order = sorted(
        ids, key=lambda i: (sp.default_sort_key(pts[i][0]), sp.default_sort_key(pts[i][1]))
    )

    def cross(i, j, k):
        p, q, r = pts[i], pts[j], pts[k]
        return sp.expand((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))

    from .exact_arithmetic import compare_exact_reals

    lo = []
    for i in order:
        while len(lo) >= 2 and compare_exact_reals(cross(lo[-2], lo[-1], i), 0) <= 0:
            lo.pop()
        lo.append(i)
    hi = []
    for i in reversed(order):
        while len(hi) >= 2 and compare_exact_reals(cross(hi[-2], hi[-1], i), 0) <= 0:
            hi.pop()
        hi.append(i)
    return lo[:-1] + hi[:-1]


def _boundary(region):
    if isinstance(region, (PolygonalSet, Polyhedron, PolyhedralShell)):
        return region
    if isinstance(region, Polygon):
        return PolygonalSet((region,))
    if isinstance(region, Polytope) and region.dimension() == 3 and region.ambient_dimension() == 3:
        facets = region.facets()
        verts = region.vertices
        center = sp.Matrix([sum(v[j] for v in verts) / len(verts) for j in range(3)])
        faces = []
        for facet in facets:
            ids = _ordered_facet(facet.vertex_indices, verts)
            # facet indices are cyclic in incidence; orient outward using Newell-like cross from first 3
            if len(ids) >= 3:
                a, b, c = (sp.Matrix(verts[i]) for i in ids[:3])
                normal = (b - a).cross(c - a)
                fc = sum((sp.Matrix(verts[i]) for i in ids), sp.zeros(3, 1)) / len(ids)
                if certified_sign(normal.dot(center - fc)) == 1:
                    ids.reverse()
            faces.append(tuple(ids))
        shell = PolyhedralShell(verts, faces)
        return Polyhedron((PolyhedralComponent(shell),))
    raise NotImplementedError(f"boundary conversion is not implemented for {type(region).__name__}")


def convert_region(
    region,
    representation="canonical",
    *,
    variables=None,
    return_result=False,
    bounds=None,
    dimension=None,
    require_verified=True,
    slices=4,
    tolerance=1e-9,
    precision=30,
    require_conforming=True,
):
    """Convert a region to a certified exact or explicitly lossy representation."""
    aliases = {"implicit": "formula", "triangulation": "simplicial", "boundary_mesh": "boundary"}
    target = aliases.get(representation, representation)
    if target == "canonical":
        value = canonicalize_region(region)
        exact = cert = True
        lossy = False
    elif target in {"formula", "semialgebraic"}:
        reg = _semialgebraic(region, variables)
        value = reg.formula if target == "formula" else reg
        exact = cert = True
        lossy = False
    elif target == "cad":
        from .cad_region import as_cad_region

        value = as_cad_region(_semialgebraic(region, variables))
        exact = cert = True
        lossy = False
    elif target == "parametric":
        if not isinstance(region, Geometry):
            raise TypeError("parametric conversion requires explicit Geometry")
        value = (
            region.bounded_parametric_cover(variables, bounds)
            if bounds is not None
            else region.intrinsic_parametric_cover(
                variables, dimension=dimension, require_verified=require_verified
            )
        )
        exact = cert = True
        lossy = False
    elif target == "simplicial":
        from .topology.semialgebraic import triangulate_region

        value = triangulate_region(region, variables)
        exact = cert = True
        lossy = False
    elif target == "mesh":
        from .cad_algorithms.meshing import triangulate_cad_region

        value = triangulate_cad_region(
            _semialgebraic(region, variables),
            slices=slices,
            tolerance=tolerance,
            precision=precision,
            require_conforming=require_conforming,
        )
        exact = False
        cert = True
        lossy = True
    elif target == "boundary":
        value = _boundary(region)
        exact = cert = True
        lossy = False
    else:
        raise ValueError("unknown representation")
    result = RegionConversion(value, target, exact, cert, lossy, type(region).__name__)
    return result if return_result else value
