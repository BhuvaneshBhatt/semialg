from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..algebraic.samples import Sample, sample_to_expr
from ..cad_algorithms.lifting.stack import CADCell
from ..context import with_computation_context
from ..domains import normalize_assumptions
from ..qe.complete import cells_to_formula
from ..simplify.formula import simplify_qe_formula
from ..structural_keys import symbol_identity_key
from ..symbol_resolution import build_symbol_table
from ..topology.incidence import cell_dimension, closures_intersect
from .cylindrical import CellSet, cad


@dataclass(frozen=True)
class CADComponent:
    """One connected component represented by selected CAD cells."""

    id: int
    variables: tuple[sp.Symbol, ...]
    cells: tuple[CADCell, ...]
    dimension: int
    formula: sp.Expr
    sample_exact: Mapping[sp.Symbol, sp.Expr]
    sample_approx: Mapping[sp.Symbol, float]

    def to_rules(self) -> Mapping[sp.Symbol, sp.Expr]:
        return self.sample_exact

    def sample_point(self) -> Mapping[sp.Symbol, sp.Expr]:
        return self.sample_exact


@dataclass(frozen=True)
class ComponentResult:
    """Result returned by :func:`component_instances`."""

    variables: tuple[sp.Symbol, ...]
    components: tuple[CADComponent, ...]
    status: str
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    @property
    def instances(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return tuple(comp.sample_exact for comp in self.components)

    @property
    def approx_instances(self) -> tuple[Mapping[sp.Symbol, float], ...]:
        return tuple(comp.sample_approx for comp in self.components)

    def __len__(self) -> int:
        return len(self.components)

    def __iter__(self):
        return iter(self.components)


@dataclass(frozen=True)
class CellGraph:
    """Undirected adjacency graph for selected CAD cells."""

    nodes: tuple[tuple[int, ...], ...]
    edges: Mapping[tuple[int, ...], frozenset[tuple[int, ...]]]

    def neighbors(self, node: tuple[int, ...]) -> frozenset[tuple[int, ...]]:
        return self.edges.get(node, frozenset())


@with_computation_context
def component_instances(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol | str],
    *,
    domain: str = "reals",
    strategy: str = "auto",
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    max_components: int | None = None,
    strict: bool = False,
    return_result: bool = False,
):
    """Return one exact sample mapping from each connected component.

    The default return is the tuple of mappings. Set ``return_result=True``
    for component objects, status, and CAD diagnostics.
    """

    result = cad(
        formula,
        variables,
        output="cells",
        domain=domain,
        strategy=strategy,
        assumptions=assumptions,
        strict=strict,
        return_result=True,
    )
    comp_result = components_from_cell_set(
        result.as_cell_set(),
        max_components=max_components,
        diagnostics={
            "cad_status": result.status,
            "strategy": strategy,
            "domain": domain,
            "assumptions": tuple(map(sp.sstr, normalize_assumptions(assumptions))),
        },
    )
    if return_result:
        return comp_result
    return comp_result.instances


def components_from_cell_set(
    cell_set: CellSet,
    *,
    max_components: int | None = None,
    diagnostics: Mapping[str, object] | None = None,
) -> ComponentResult:
    """Build connected components from an existing CAD cell set."""

    graph = cell_adjacency_graph(cell_set)
    index_to_cell = {cell.index: cell for cell in cell_set.cells}
    parts = _connected_parts(graph)
    if max_components is not None:
        parts = parts[:max_components]
    comps: list[CADComponent] = []
    for comp_id, part in enumerate(parts):
        cells = tuple(index_to_cell[index] for index in sorted(part))
        comps.append(_make_component(comp_id, cells, cell_set))
    merged_diag = dict(diagnostics or {})
    merged_diag.update(
        {
            "cell_count": len(cell_set.cells),
            "component_count": len(comps),
            "graph_edge_count": sum(len(v) for v in graph.edges.values()) // 2,
        }
    )
    return ComponentResult(
        variables=cell_set.variables,
        components=tuple(comps),
        status="complete",
        diagnostics=merged_diag,
    )


def cell_adjacency_graph(cell_set: CellSet) -> CellGraph:
    """Construct the closure-connected graph induced by selected CAD cells."""

    cells = tuple(cell_set.cells)
    final_level = len(cell_set.variables)
    candidates = tuple(cell_set.cells_by_level.get(final_level, ()))
    edges: dict[tuple[int, ...], set[tuple[int, ...]]] = {cell.index: set() for cell in cells}
    for pos, left in enumerate(cells):
        for right in cells[pos + 1 :]:
            if len(cell_set.variables) == 1:
                adjacent = abs(left.stack_position - right.stack_position) <= 1
            else:
                adjacent = closures_intersect(
                    left, right, cell_set.variables, cell_set.cells_by_level, candidates
                )
            if adjacent:
                edges[left.index].add(right.index)
                edges[right.index].add(left.index)
    frozen = {node: frozenset(neighbors) for node, neighbors in edges.items()}
    return CellGraph(nodes=tuple(edges), edges=frozen)


def _connected_parts(graph: CellGraph) -> list[frozenset[tuple[int, ...]]]:
    unseen = set(graph.nodes)
    parts: list[frozenset[tuple[int, ...]]] = []
    while unseen:
        start = min(unseen)
        stack = [start]
        seen = {start}
        unseen.remove(start)
        while stack:
            node = stack.pop()
            for neighbor in graph.neighbors(node):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    seen.add(neighbor)
                    stack.append(neighbor)
        parts.append(frozenset(seen))
    return parts


def _make_component(comp_id: int, cells: tuple[CADCell, ...], cell_set: CellSet) -> CADComponent:
    dim = max((cell_dimension(cell, cell_set.cells_by_level) for cell in cells), default=0)
    sample_cell = _best_sample_cell(cells, cell_set)
    exact = _sample_map(cell_set.variables, sample_cell.sample if sample_cell is not None else ())
    approx = {var: float(sp.N(value)) for var, value in exact.items()}
    formula = (
        simplify_qe_formula(
            cells_to_formula(cells, cell_set.variables, cell_set.cells_by_level),
            implication_minimize=False,
        )
        if cells
        else sp.false
    )
    return CADComponent(
        id=comp_id,
        variables=cell_set.variables,
        cells=cells,
        dimension=dim,
        formula=formula,
        sample_exact=exact,
        sample_approx=approx,
    )


def _best_sample_cell(cells: Sequence[CADCell], cell_set: CellSet) -> CADCell | None:
    if not cells:
        return None
    return max(
        cells, key=lambda cell: (cell_dimension(cell, cell_set.cells_by_level), -sum(cell.index))
    )


def _sample_map(
    variables: Sequence[sp.Symbol], samples: Sequence[Sample]
) -> Mapping[sp.Symbol, sp.Expr]:
    return {var: sample_to_expr(sample) for var, sample in zip(variables, samples, strict=True)}


@with_computation_context
def component_instances_text(
    text: str,
    *,
    variables: Sequence[sp.Symbol | str] | None = None,
    symbols: Mapping[str, sp.Symbol] | None = None,
    domain: str = "reals",
    strategy: str = "auto",
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    max_components: int | None = None,
    strict: bool = False,
    return_result: bool = False,
):
    """Text wrapper for :func:`component_instances`."""

    from ..formula import parse_formula_text
    from .cylindrical import _normalize_variables

    local_symbols = build_symbol_table(symbols, variables or ())
    expr, _ = parse_formula_text(text, symbols=local_symbols)
    var_tuple = (
        _normalize_variables(variables)
        if variables is not None
        else tuple(sorted(expr.free_symbols, key=symbol_identity_key))
    )
    return component_instances(
        expr,
        var_tuple,
        domain=domain,
        strategy=strategy,
        assumptions=assumptions,
        max_components=max_components,
        strict=strict,
        return_result=return_result,
    )


__all__ = [
    "CADComponent",
    "CellGraph",
    "ComponentResult",
    "cell_adjacency_graph",
    "component_instances",
    "component_instances_text",
    "components_from_cell_set",
]
