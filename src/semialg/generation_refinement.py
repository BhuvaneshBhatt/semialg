"""Exact/deterministic polyhedral generation and triangular-face refinement."""

from __future__ import annotations

import random

import sympy as sp

from ._zero_testing import certified_zero


def random_polytope(dimension=3, *, point_count=None, coordinate_bound=10, seed=None):
    """Generate a reproducible exact full-dimensional random lattice polytope."""
    if dimension < 1:
        raise ValueError("dimension must be positive")
    if coordinate_bound < 1:
        raise ValueError("coordinate_bound must be positive")
    count = point_count if point_count is not None else max(dimension + 1, 2 * dimension + 2)
    if count < dimension + 1:
        raise ValueError("point_count must be at least dimension + 1")
    rng = random.Random(seed)
    for _ in range(128):
        pts = tuple(
            tuple(
                sp.Integer(rng.randint(-coordinate_bound, coordinate_bound))
                for _ in range(dimension)
            )
            for _ in range(count)
        )
        if len(set(pts)) < dimension + 1:
            continue
        anchor = sp.Matrix(pts[0])
        if sp.Matrix.hstack(*(sp.Matrix(p) - anchor for p in pts[1:])).rank() != dimension:
            continue
        from .convex_hull import convex_hull

        return convex_hull(pts)
    raise RuntimeError("could not generate a full-dimensional random polytope")


def random_polygon(*, vertex_count=6, coordinate_bound=10, seed=None):
    """Generate a reproducible exact convex lattice polygon."""
    if vertex_count < 3:
        raise ValueError("vertex_count must be at least three")
    # Oversample because interior lattice points disappear under the hull.
    rng = random.Random(seed)
    from .convex_hull import convex_hull
    from .standard_regions import Polygon, Polytope

    for _ in range(128):
        pts = tuple(
            (
                sp.Integer(rng.randint(-coordinate_bound, coordinate_bound)),
                sp.Integer(rng.randint(-coordinate_bound, coordinate_bound)),
            )
            for _ in range(max(vertex_count * 3, 12))
        )
        hull = convex_hull(pts, canonical=False)
        if isinstance(hull, Polytope) and len(hull.vertices) >= 3:
            incidence = hull.incidence()
            edges = [facet.vertex_indices for facet in incidence.facets]
            neighbors = {i: [] for i in range(len(hull.vertices))}
            for i, j in edges:
                neighbors[i].append(j)
                neighbors[j].append(i)
            if all(len(value) == 2 for value in neighbors.values()):
                start = min(neighbors)
                ordered_indices = [start]
                previous = None
                current = start
                while True:
                    choices = sorted(neighbors[current])
                    nxt = choices[0] if choices[0] != previous else choices[1]
                    if nxt == start:
                        break
                    ordered_indices.append(nxt)
                    previous, current = current, nxt
                ordered = [hull.vertices[i] for i in ordered_indices]
                if len(ordered) >= vertex_count:
                    return Polygon(ordered)
    raise RuntimeError("could not generate a random polygon")


def _triangular_shell(source):
    from .boundary_topology import PolyhedralShell
    from .standard_regions import Polytope

    if isinstance(source, PolyhedralShell):
        shell = source
    elif isinstance(source, Polytope):
        from .region_conversion import convert_region

        boundary = convert_region(source, "boundary")
        shell = boundary.components[0].outer_shell if hasattr(boundary, "components") else boundary
    else:
        raise TypeError("source must be a Polytope or PolyhedralShell")
    if any(len(face) != 3 for face in shell.faces):
        raise ValueError("triangular-face refinement requires a triangular shell")
    return shell


def subdivide_triangular_faces(source, *, levels=1):
    """Subdivide every triangular face into four triangles using shared exact midpoints."""
    from .boundary_topology import PolyhedralShell

    if levels < 0:
        raise ValueError("levels must be nonnegative")
    shell = _triangular_shell(source)
    for _ in range(levels):
        vertices = list(shell.vertices)
        midpoint_index = {}

        def midpoint(i, j, *, _vertices=vertices, _indices=midpoint_index):
            edge = tuple(sorted((i, j)))
            if edge not in _indices:
                p, q = _vertices[i], _vertices[j]
                m = tuple(sp.simplify((a + b) / 2) for a, b in zip(p, q, strict=True))
                _indices[edge] = len(_vertices)
                _vertices.append(m)
            return _indices[edge]

        faces = []
        for a, b, c in shell.faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            faces.extend(((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)))
        shell = PolyhedralShell(vertices, faces)
    return shell


def geodesic_refinement(source, *, levels=1, center=(0, 0, 0), radius=None):
    """Refine a triangular shell and project its vertices exactly to a sphere."""
    from .boundary_topology import PolyhedralShell

    shell = subdivide_triangular_faces(source, levels=levels)
    c = tuple(sp.sympify(x) for x in center)
    if len(c) != 3:
        raise ValueError("center must be three-dimensional")
    if radius is None:
        p0 = shell.vertices[0]
        radius = sp.sqrt(sum((a - b) ** 2 for a, b in zip(p0, c, strict=True)))
    r = sp.sympify(radius)
    if certified_zero(r) is True:
        raise ValueError("radius must be nonzero")
    projected = []
    for p in shell.vertices:
        delta = tuple(sp.simplify(a - b) for a, b in zip(p, c, strict=True))
        norm = sp.sqrt(sum(x * x for x in delta))
        if certified_zero(norm) is True:
            raise ValueError("cannot project the sphere center")
        projected.append(
            tuple(sp.simplify(a + r * x / norm) for a, x in zip(c, delta, strict=True))
        )
    return PolyhedralShell(projected, shell.faces)


__all__ = ["random_polygon", "random_polytope", "subdivide_triangular_faces", "geodesic_refinement"]
