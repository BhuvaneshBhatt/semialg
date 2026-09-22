"""Conforming exact tetrahedralization of mixed 3D convex cells."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .polytope_decomposition import triangulate_polytope
from .standard_regions import Polytope


@dataclass(frozen=True)
class MixedCellTetrahedralization:
    """Certified tetrahedralization of one convex three-dimensional cell."""

    vertex_ids: tuple[int, ...]
    tetrahedra: tuple[tuple[int, int, int, int], ...]
    certified: bool = True
    conforming: bool = True


@dataclass(frozen=True)
class MixedMeshTetrahedralization:
    """Conforming tetrahedralization of a mixed convex-cell mesh."""

    vertices: tuple[tuple[sp.Expr, ...], ...]
    tetrahedra: tuple[tuple[int, int, int, int], ...]
    source_cell_ranges: tuple[tuple[int, int], ...]
    certified: bool = True
    conforming: bool = True


def tetrahedralize_cell(
    vertices: Sequence[Sequence[object]], *, vertex_ids: Sequence[int] | None = None
) -> MixedCellTetrahedralization:
    """Tetrahedralize one convex 3-cell using a global vertex ordering."""
    pts = tuple(tuple(sp.sympify(x) for x in p) for p in vertices)
    ids = tuple(range(len(pts))) if vertex_ids is None else tuple(map(int, vertex_ids))
    if len(ids) != len(pts) or len(set(ids)) != len(ids):
        raise ValueError("vertex_ids must be distinct and match vertices")
    order = sorted(range(len(ids)), key=lambda i: ids[i])
    ordered = tuple(pts[i] for i in order)
    dec = triangulate_polytope(Polytope(ordered), strategy="pulling")
    local_to_global = tuple(ids[i] for i in order)
    tets = tuple(tuple(local_to_global[i] for i in s) for s in dec.simplices)
    if any(len(t) != 4 for t in tets):
        raise ValueError("cell must be three-dimensional")
    return MixedCellTetrahedralization(ids, tets)


def tetrahedralize_cells(
    vertices: Sequence[Sequence[object]], cells: Sequence[Sequence[int]]
) -> MixedMeshTetrahedralization:
    """Tetrahedralize cells conformingly using shared global vertex identifiers."""
    pts = tuple(tuple(sp.sympify(x) for x in p) for p in vertices)
    all_tets = []
    ranges = []
    normalized_cells = []
    for cell in cells:
        ids = tuple(map(int, cell))
        if len(set(ids)) != len(ids):
            raise ValueError("a cell cannot repeat a global vertex identifier")
        if any(i < 0 or i >= len(pts) for i in ids):
            raise IndexError("cell vertex index is outside vertex array")
        start = len(all_tets)
        result = tetrahedralize_cell(tuple(pts[i] for i in ids), vertex_ids=ids)
        all_tets.extend(result.tetrahedra)
        ranges.append((start, len(all_tets)))
        normalized_cells.append(ids)
    # Replay the conformity contract on every shared polygonal face candidate:
    # the tetrahedral boundary induced on the common global vertices must agree.
    import itertools

    def induced_triangles(start, end, common):
        counts = {}
        for tet in all_tets[start:end]:
            for tri in itertools.combinations(tet, 3):
                key = frozenset(tri)
                counts[key] = counts.get(key, 0) + 1
        return {tri for tri, count in counts.items() if count == 1 and tri <= common}

    for i in range(len(normalized_cells)):
        for j in range(i + 1, len(normalized_cells)):
            common = frozenset(normalized_cells[i]) & frozenset(normalized_cells[j])
            if len(common) < 3:
                continue
            left = induced_triangles(*ranges[i], common)
            right = induced_triangles(*ranges[j], common)
            if left != right:
                raise ValueError(
                    "adjacent cells induce incompatible triangulations on a shared face"
                )
    return MixedMeshTetrahedralization(pts, tuple(all_tets), tuple(ranges))
