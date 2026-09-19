"""Numerical simplicial meshing for structured CAD cells.

This module contains numerical geometry only; exact reusable CAD-region algebra
remains in :mod:`semialg.cad_region`.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from .numerical_boundaries import _numeric_real_roots


@dataclass(frozen=True)
class CADMesh:
    """Numerical simplicial mesh carrying exact CAD-cell provenance."""

    coordinates: tuple[tuple[float, ...], ...]
    simplices: tuple[tuple[int, ...], ...]
    cell_index: tuple[int, ...] = ()
    simplex_cells: tuple[tuple[int, ...], ...] = ()
    vertex_cells: tuple[tuple[tuple[int, ...], ...], ...] = ()
    conforming: bool | None = None

    @property
    def dimension(self) -> int:
        return len(self.simplices[0]) - 1 if self.simplices else 0

    @property
    def shared_vertex_count(self) -> int:
        return sum(1 for owners in self.vertex_cells if len(owners) > 1)

    @property
    def source_cells(self) -> tuple[tuple[int, ...], ...]:
        return tuple(sorted(set(self.simplex_cells)))


def _float_real(expr: sp.Expr, subs: Mapping[sp.Symbol, object] | None = None) -> float:
    val = sp.N(expr.subs(subs or {}), 30)
    if val.is_real is False:
        raise ValueError(f"non-real CAD bound value: {expr}")
    out = float(val)
    if not math.isfinite(out):
        raise ValueError(f"non-finite CAD bound value: {expr}")
    return out


def _numeric_bound_value(
    level, side: str, subs: Mapping[sp.Symbol, float], *, precision: int = 30
) -> float:
    """Evaluate one finite typed CAD bound, preserving its root index."""
    from .bounds import AlgebraicRootFunction

    bound = level.lower_bound if side == "lower" else level.upper_bound
    expr = level.lower if side == "lower" else level.upper
    if isinstance(bound, AlgebraicRootFunction):
        roots = _numeric_real_roots(
            bound.polynomial, bound.fiber_variable, subs, precision=precision
        )
        if bound.root_index >= len(roots):
            raise ValueError(
                f"CAD {side} boundary root {bound.root_index + 1} disappeared during meshing"
            )
        return roots[bound.root_index]
    return _float_real(expr, subs)


def _unit_grid_indices(shape: Sequence[int]) -> Iterable[tuple[int, ...]]:
    import itertools

    return itertools.product(*(range(n + 1) for n in shape))


def _grid_vertex_id(index: Sequence[int], shape: Sequence[int]) -> int:
    out = 0
    stride = 1
    for i, n in zip(reversed(index), reversed(shape), strict=True):
        out += i * stride
        stride *= n + 1
    return out


def _freudenthal_simplices(shape: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """Triangulate a rectangular integer grid by the Freudenthal scheme."""
    import itertools

    d = len(shape)
    simplices: list[tuple[int, ...]] = []
    for base in itertools.product(*(range(n) for n in shape)):
        for perm in itertools.permutations(range(d)):
            vertex = list(base)
            ids = [_grid_vertex_id(vertex, shape)]
            for axis in perm:
                vertex = list(vertex)
                vertex[axis] += 1
                ids.append(_grid_vertex_id(vertex, shape))
            simplices.append(tuple(ids))
    return tuple(simplices)


def _triangulate_planar_strip(cell, *, slices: int, inset: float) -> CADMesh:
    """Preserve the original strip triangulation and its compact vertex count."""
    xlevel, ylevel = cell.levels
    if xlevel.lower == -sp.oo or xlevel.upper == sp.oo:
        raise ValueError("CAD cell must have finite x bounds")
    x0, x1 = _float_real(xlevel.lower), _float_real(xlevel.upper)
    if not x0 < x1:
        raise ValueError("CAD cell has degenerate x interval")
    eps = min(abs(x1 - x0) * inset, abs(x1 - x0) / (10 * slices))
    xs = [x0 + eps + (x1 - x0 - 2 * eps) * i / slices for i in range(slices + 1)]
    coords: list[tuple[float, float]] = []
    xvar = cell.variables[0]
    for xv in xs:
        lo = _numeric_bound_value(ylevel, "lower", {xvar: xv})
        hi = _numeric_bound_value(ylevel, "upper", {xvar: xv})
        if not lo < hi:
            raise ValueError("CAD cell vertical bounds cross or collapse numerically")
        coords.extend(((xv, lo), (xv, hi)))
    tris: list[tuple[int, int, int]] = []
    for i in range(slices):
        a, b, c, d = 2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 3
        tris.extend(((a, c, d), (a, d, b)))
    idx = tuple(cell.index)
    return CADMesh(tuple(coords), tuple(tris), idx, tuple(idx for _ in tris))


def triangulate_cad_cell(
    cell: object,
    *,
    slices: int | Sequence[int] = 12,
    inset: float = 1e-8,
    precision: int = 30,
) -> CADMesh:
    """Numerically simplicially mesh a bounded full-dimensional CAD cell.

    Two-dimensional cells use a strip triangulation. In dimensions >=3 a
    tensor grid on the unit cube is mapped
    recursively through the cylindrical bounds and each grid cube is split by
    the Freudenthal triangulation.  The mapping evaluates typed algebraic root
    functions by their CAD root indices rather than re-solving for a branch.
    """
    from .structured_cells import StructuredCADCell

    if not isinstance(cell, StructuredCADCell):
        raise TypeError("triangulate_cad_cell expects a StructuredCADCell")
    d = len(cell.variables)
    if d < 2 or not cell.is_full_dimensional:
        raise ValueError("triangulation requires a full-dimensional CAD cell of dimension >= 2")
    if not cell.bounded:
        raise ValueError("CAD cell must have finite bounds at every coordinate level")
    if isinstance(slices, int):
        if slices < 1:
            raise ValueError("slices must be positive")
        if d == 2:
            if slices < 2:
                raise ValueError("slices must be at least 2 in dimension two")
            return _triangulate_planar_strip(cell, slices=slices, inset=inset)
        shape = (slices,) * d
    else:
        shape = tuple(int(n) for n in slices)
        if len(shape) != d or any(n < 1 for n in shape):
            raise ValueError("slices sequence must have one positive entry per dimension")
        if d == 2 and shape[1] == 1:
            return _triangulate_planar_strip(cell, slices=shape[0], inset=inset)

    coords: list[tuple[float, ...]] = []
    for grid_index in _unit_grid_indices(shape):
        physical: list[float] = []
        subs: dict[sp.Symbol, float] = {}
        for axis, (level, n, gi) in enumerate(zip(cell.levels, shape, grid_index, strict=True)):
            lo = _numeric_bound_value(level, "lower", subs, precision=precision)
            hi = _numeric_bound_value(level, "upper", subs, precision=precision)
            if not lo < hi:
                raise ValueError(f"CAD bounds cross or collapse numerically at level {axis + 1}")
            # Meshing a cell represents its closure.  inset>0 is available when
            # callers explicitly want an interior-only numerical approximation.
            t = gi / n
            if inset > 0:
                eps = min(float(inset), 0.25 / n)
                t = eps + (1.0 - 2.0 * eps) * t
            value = lo + (hi - lo) * t
            physical.append(value)
            subs[cell.variables[axis]] = value
        coords.append(tuple(physical))
    simplices = _freudenthal_simplices(shape)
    idx = tuple(cell.index)
    return CADMesh(tuple(coords), simplices, idx, tuple(idx for _ in simplices))


def _merge_vertex(
    point: tuple[float, ...],
    coordinates: list[tuple[float, ...]],
    buckets: dict[tuple[int, ...], list[int]],
    *,
    tolerance: float,
) -> int:
    if tolerance <= 0:
        try:
            return coordinates.index(point)
        except ValueError:
            coordinates.append(point)
            return len(coordinates) - 1
    key = tuple(round(x / tolerance) for x in point)
    # Inspect neighboring quantization buckets so points straddling a rounding
    # boundary still merge.
    import itertools

    for delta in itertools.product((-1, 0, 1), repeat=len(point)):
        near = tuple(k + dk for k, dk in zip(key, delta, strict=True))
        for idx in buckets.get(near, ()):
            old = coordinates[idx]
            if max(abs(a - b) for a, b in zip(old, point, strict=True)) <= tolerance:
                return idx
    idx = len(coordinates)
    coordinates.append(point)
    buckets.setdefault(key, []).append(idx)
    return idx


def _mesh_is_conforming(
    simplices: Sequence[tuple[int, ...]],
    simplex_cells: Sequence[tuple[int, ...]],
    vertex_cells: Sequence[Sequence[tuple[int, ...]]],
) -> bool:
    """Check that adjacent source cells induce identical shared simplex facets."""
    from itertools import combinations

    if not simplices:
        return True
    facets_by_owner: dict[tuple[int, ...], set[tuple[int, ...]]] = {}
    for simplex, owner in zip(simplices, simplex_cells, strict=True):
        for drop in range(len(simplex)):
            facet = tuple(sorted(v for i, v in enumerate(simplex) if i != drop))
            facets_by_owner.setdefault(owner, set()).add(facet)
    owners = tuple(sorted(facets_by_owner))
    for a, b in combinations(owners, 2):
        shared_vertices = {
            i for i, provenance in enumerate(vertex_cells) if a in provenance and b in provenance
        }
        if not shared_vertices:
            continue
        af = {f for f in facets_by_owner[a] if set(f) <= shared_vertices}
        bf = {f for f in facets_by_owner[b] if set(f) <= shared_vertices}
        if af != bf:
            return False
    return True


def triangulate_cad_cells(
    cells: Iterable[object],
    *,
    slices: int | Sequence[int] = 4,
    tolerance: float = 1e-9,
    precision: int = 30,
    require_conforming: bool = True,
) -> CADMesh:
    """Mesh multiple adjacent CAD cells with one shared vertex table.

    Individual cells are meshed on their closures (``inset=0``).  Geometrically
    coincident vertices are canonicalized, so neighboring CAD cells reference
    the same vertex indices on sampled common faces.  ``simplex_cells`` and
    ``vertex_cells`` retain exact CAD-index provenance.
    """
    from .structured_cells import StructuredCADCell

    items = tuple(cells)
    if not items:
        return CADMesh((), ())
    if any(not isinstance(cell, StructuredCADCell) for cell in items):
        raise TypeError("triangulate_cad_cells expects StructuredCADCell objects")
    dimensions = {len(cell.variables) for cell in items}
    if len(dimensions) != 1:
        raise ValueError("all CAD cells must have the same ambient dimension")
    coordinates: list[tuple[float, ...]] = []
    buckets: dict[tuple[int, ...], list[int]] = {}
    simplices: list[tuple[int, ...]] = []
    simplex_cells: list[tuple[int, ...]] = []
    vertex_provenance: list[set[tuple[int, ...]]] = []
    for cell in items:
        local = triangulate_cad_cell(cell, slices=slices, inset=0.0, precision=precision)
        remap: list[int] = []
        for point in local.coordinates:
            before = len(coordinates)
            global_idx = _merge_vertex(point, coordinates, buckets, tolerance=tolerance)
            if len(coordinates) > before:
                vertex_provenance.append(set())
            vertex_provenance[global_idx].add(tuple(cell.index))
            remap.append(global_idx)
        for simplex in local.simplices:
            mapped = tuple(remap[i] for i in simplex)
            if len(set(mapped)) == len(mapped):
                simplices.append(mapped)
                simplex_cells.append(tuple(cell.index))
    vertex_cells = tuple(tuple(sorted(cells_)) for cells_ in vertex_provenance)
    conforming = _mesh_is_conforming(simplices, simplex_cells, vertex_cells)
    if require_conforming and not conforming:
        raise ValueError(
            "adjacent CAD cells produced incompatible shared-face triangulations; "
            "use a common slice grid or set require_conforming=False to inspect the raw mesh"
        )
    return CADMesh(
        tuple(coordinates),
        tuple(simplices),
        (),
        tuple(simplex_cells),
        vertex_cells,
        conforming,
    )


def triangulate_cad_region(
    region: object,
    *,
    slices: int | Sequence[int] = 4,
    tolerance: float = 1e-9,
    precision: int = 30,
    require_conforming: bool = True,
) -> CADMesh:
    """Adjacency-aware mesh of all bounded selected full-dimensional CAD cells."""
    from ..cad_region import as_cad_region
    from .structured_cells import extract_structured_cad_cells

    reg = as_cad_region(region)
    decomposition = extract_structured_cad_cells(reg.result)
    cells = tuple(cell for cell in decomposition.full_dimensional_cells if cell.bounded)
    return triangulate_cad_cells(
        cells,
        slices=slices,
        tolerance=tolerance,
        precision=precision,
        require_conforming=require_conforming,
    )


__all__ = ["CADMesh", "triangulate_cad_cell", "triangulate_cad_cells", "triangulate_cad_region"]
