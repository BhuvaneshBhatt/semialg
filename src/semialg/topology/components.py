from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

from ..cad.decomposition import CompleteCAD
from ..cad.lifting.stack import CADCell
from .incidence import closures_intersect, final_cells


def component_cell_graph(
    selected: Iterable[CADCell],
    cad: CompleteCAD,
    variables: Sequence,
) -> dict[tuple[int, ...], frozenset[tuple[int, ...]]]:
    """Build the closure-connected adjacency graph induced by selected cells."""

    cells = tuple(sorted(selected, key=lambda cell: cell.index))
    candidate_boundary = final_cells(cad)
    edges: dict[tuple[int, ...], set[tuple[int, ...]]] = {cell.index: set() for cell in cells}
    for pos, left in enumerate(cells):
        for right in cells[pos + 1 :]:
            if closures_intersect(left, right, variables, cad.cells_by_level, candidate_boundary):
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
