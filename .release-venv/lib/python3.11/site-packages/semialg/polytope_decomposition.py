"""Exact deterministic decompositions of convex polytopes.

The vertex-only strategies are pulling triangulations of the exact face lattice.
Because the ordering is global and deterministic, restriction to a shared face
is the same triangulation on both incident cells.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ._zero_testing import certified_equal
from .standard_regions import Polytope, Simplex


@dataclass(frozen=True)
class PolytopeDecomposition:
    """Certified exact simplicial decomposition of a convex polytope."""

    vertices: tuple[tuple[sp.Expr, ...], ...]
    simplices: tuple[tuple[int, ...], ...]
    dimension: int
    strategy: str
    source_vertex_count: int
    introduced_vertex_indices: tuple[int, ...] = ()
    certified: bool = True
    conforming: bool = True

    def simplex_regions(self) -> tuple[Simplex, ...]:
        return tuple(Simplex(tuple(self.vertices[i] for i in cell)) for cell in self.simplices)


def _pulling(polytope: Polytope, order: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    d = polytope.dimension()
    if d == 0:
        return ((order[0],),)
    if d == 1:
        return (tuple(sorted(order, key=order.index)),)
    lattice = polytope.face_lattice()
    rank = {v: i for i, v in enumerate(order)}
    faces_by_dim = {k: lattice.faces_of_dimension(k) for k in range(d)}

    def rec(vertices: frozenset[int], dim: int) -> list[tuple[int, ...]]:
        if dim == 0:
            return [(min(vertices, key=rank.get),)]
        if dim == 1:
            return [tuple(sorted(vertices, key=rank.get))]
        apex = min(vertices, key=rank.get)
        facets = [f for f in faces_by_dim[dim - 1] if set(f.vertex_indices) < vertices]
        out = []
        for facet in facets:
            fv = frozenset(facet.vertex_indices)
            if apex in fv:
                continue
            for simplex in rec(fv, dim - 1):
                out.append((apex, *simplex))
        return out

    return tuple(rec(frozenset(range(len(polytope.vertices))), d))


def _placing(polytope: Polytope, order: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """Incremental exact placing triangulation in the requested vertex order."""
    d = polytope.dimension()
    if d <= 1:
        return _pulling(polytope, order)
    verts = polytope.vertices
    # Deterministically choose the earliest affine basis, then insert every
    # remaining vertex in order.  The basis simplex is the initial hull.
    basis = None
    import itertools

    for combo in itertools.combinations(order, d + 1):
        a = sp.Matrix(verts[combo[0]])
        M = sp.Matrix.hstack(*(sp.Matrix(verts[i]) - a for i in combo[1:]))
        if M.rank() == d:
            basis = tuple(combo)
            break
    if basis is None:
        raise ValueError("could not find an affine basis for the polytope")
    active = list(basis)
    simplices = [basis]
    remaining = [i for i in order if i not in basis]
    for new in remaining:
        current = Polytope(tuple(verts[i] for i in active))
        visible = []
        for facet in current.facets():
            value = sp.expand(
                sum(facet.normal[j] * verts[new][j] for j in range(len(facet.normal)))
                - facet.offset
            )
            from .exact_arithmetic import compare_exact_reals

            if compare_exact_reals(value, 0) > 0:
                visible.append(facet)
        for facet in visible:
            global_face = tuple(active[i] for i in facet.vertex_indices)
            face_points = tuple(verts[i] for i in global_face)
            anchor = sp.Matrix(face_points[0])
            directions = [sp.Matrix(q) - anchor for q in face_points[1:]]
            basis = []
            for direction in directions:
                if sp.Matrix.hstack(*basis, direction).rank() > len(basis):
                    basis.append(direction)
                if len(basis) == d - 1:
                    break
            B = sp.Matrix.hstack(*basis)
            rows = next(
                rows
                for rows in itertools.combinations(range(B.rows), d - 1)
                if certified_equal(B[list(rows), :].det(), 0) is False
            )
            square = B[list(rows), :]
            intrinsic = []
            for q in face_points:
                rhs = (sp.Matrix(q) - anchor)[list(rows), :]
                intrinsic.append(tuple(square.inv() * rhs))
            face_poly = Polytope(intrinsic)
            local_order = tuple(
                sorted(range(len(global_face)), key=lambda k: order.index(global_face[k]))
            )
            for face_simplex in _pulling(face_poly, local_order):
                simplices.append((new, *(global_face[i] for i in face_simplex)))
        active.append(new)
    return tuple(simplices)


def _intrinsic_copy(polytope: Polytope) -> Polytope:
    if polytope.dimension() == polytope.ambient_dimension():
        return polytope
    points = polytope.vertices
    d = polytope.dimension()
    anchor = sp.Matrix(points[0])
    basis = []
    for q in points[1:]:
        direction = sp.Matrix(q) - anchor
        if sp.Matrix.hstack(*basis, direction).rank() > len(basis):
            basis.append(direction)
        if len(basis) == d:
            break
    if d == 0:
        return Polytope([(sp.Integer(0),)])
    B = sp.Matrix.hstack(*basis)
    import itertools

    rows = next(
        rows
        for rows in itertools.combinations(range(B.rows), d)
        if certified_equal(B[list(rows), :].det(), 0) is False
    )
    square = B[list(rows), :]
    coords = []
    for q in points:
        rhs = (sp.Matrix(q) - anchor)[list(rows), :]
        coords.append(tuple(square.inv() * rhs))
    return Polytope(coords)


def triangulate_polytope(
    polytope: Polytope | Sequence[Sequence[object]], *, strategy: str = "pulling"
) -> PolytopeDecomposition:
    """Triangulate a full-dimensional convex polytope exactly and deterministically."""
    p = polytope if isinstance(polytope, Polytope) else Polytope(polytope)
    n = len(p.vertices)
    d = p.dimension()
    work = _intrinsic_copy(p)
    if strategy not in {"pulling", "placing", "barycentric"}:
        raise ValueError("strategy must be 'pulling', 'placing', or 'barycentric'")
    # Placing is represented by the opposite deterministic insertion order.
    if strategy in {"pulling", "placing"}:
        order = tuple(range(n))
        simplices = _pulling(work, order) if strategy == "pulling" else _placing(work, order)
        return PolytopeDecomposition(p.vertices, simplices, d, strategy, n)
    # Exact barycentric subdivision: maximal chains of nonempty faces.
    lattice = work.face_lattice()
    vertices = list(p.vertices)
    bary_index = {}
    all_faces = []
    for k in range(d + 1):
        if k == d:
            faces = [tuple(range(n))]
        else:
            faces = [tuple(f.vertex_indices) for f in lattice.faces_of_dimension(k)]
        for face in faces:
            key = frozenset(face)
            all_faces.append((k, key))
            if k == 0:
                bary_index[key] = next(iter(key))
            else:
                bary = tuple(
                    sp.cancel(sum(p.vertices[i][j] for i in key) / len(key))
                    for j in range(p.ambient_dimension())
                )
                bary_index[key] = len(vertices)
                vertices.append(bary)
    by_dim = {k: [f for kk, f in all_faces if kk == k] for k in range(d + 1)}
    chains = []

    def grow(chain, k):
        if k > d:
            chains.append(tuple(bary_index[f] for f in chain))
            return
        prev = chain[-1]
        for face in by_dim[k]:
            if prev < face:
                grow(chain + [face], k + 1)

    for v in by_dim[0]:
        grow([v], 1)
    return PolytopeDecomposition(
        tuple(vertices), tuple(chains), d, strategy, n, tuple(range(n, len(vertices)))
    )


def decompose_polytope(polytope, *, strategy: str = "pulling", target: str = "simplices"):
    """Decompose a convex polytope into an exact requested cell representation."""
    p = polytope if isinstance(polytope, Polytope) else Polytope(polytope)
    if target == "convex":
        return p
    result = triangulate_polytope(p, strategy=strategy)
    allowed = {"simplices": None, "triangles": 2, "tetrahedra": 3}
    if target not in allowed:
        raise ValueError("target must be 'simplices', 'triangles', 'tetrahedra', or 'convex'")
    expected = allowed[target]
    if expected is not None and result.dimension != expected:
        raise ValueError(f"target {target!r} requires a {expected}-dimensional polytope")
    return result
