from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import sympy as sp

from ..cad_algorithms.decomposition import CompleteCAD
from ..cad_algorithms.lifting.stack import CADCell
from ..reconstruct.cylindrical import path_condition
from .incidence import _truth_condition_at_assignments, cell_dimension, cell_sample_subs


def component_cell_graph(
    selected: Iterable[CADCell],
    cad: CompleteCAD,
    variables: Sequence,
    *,
    closed_paths: Mapping[tuple[int, ...], sp.Expr] | None = None,
    sample_assignments: Mapping[tuple[int, ...], Mapping[sp.Symbol, sp.Expr]] | None = None,
) -> dict[tuple[int, ...], frozenset[tuple[int, ...]]]:
    """Connect selected cells through exact closure intersections.

    Closure membership is computed once per ordered source/target pair and
    reused for every cell-pair intersection.  This avoids rebuilding the same
    closed cylindrical path and re-evaluating it for each candidate common
    boundary cell.
    """

    cells = tuple(sorted(selected, key=lambda cell: cell.index))
    edges: dict[tuple[int, ...], set[tuple[int, ...]]] = {cell.index: set() for cell in cells}
    dimensions = {cell.index: cell_dimension(cell, cad.cells_by_level) for cell in cells}
    sample_assignments = dict(
        sample_assignments or {cell.index: cell_sample_subs(cell, variables) for cell in cells}
    )
    closure_members: dict[tuple[int, ...], frozenset[tuple[int, ...]]] = {}
    for source in cells:
        source_dim = dimensions[source.index]
        condition = (closed_paths or {}).get(source.index)
        if condition is None:
            condition = path_condition(source, variables, cad.cells_by_level, closed=True)
        members = {source.index}
        for target in cells:
            if target.index == source.index or dimensions[target.index] > source_dim:
                continue
            try:
                incident = _truth_condition_at_assignments(
                    condition, sample_assignments[target.index]
                )
            except (NotImplementedError, ValueError, sp.PolynomialError):
                incident = False
            if incident:
                members.add(target.index)
        closure_members[source.index] = frozenset(members)

    for pos, left in enumerate(cells):
        left_closure = closure_members[left.index]
        for right in cells[pos + 1 :]:
            if left_closure.intersection(closure_members[right.index]):
                edges[left.index].add(right.index)
                edges[right.index].add(left.index)
    return {node: frozenset(neighbors) for node, neighbors in edges.items()}


def connected_cell_components(
    graph: Mapping[tuple[int, ...], frozenset[tuple[int, ...]]],
) -> tuple[frozenset[tuple[int, ...]], ...]:
    unseen = set(graph)
    parts: list[frozenset[tuple[int, ...]]] = []
    while unseen:
        start = min(unseen)
        stack = [start]
        seen = {start}
        unseen.remove(start)
        while stack:
            node = stack.pop()
            for neighbor in graph.get(node, frozenset()):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    seen.add(neighbor)
                    stack.append(neighbor)
        parts.append(frozenset(seen))
    return tuple(parts)


__all__ = ["component_cell_graph", "connected_cell_components"]
