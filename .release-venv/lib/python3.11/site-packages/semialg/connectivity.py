from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from .cad_algorithms.cells import (
    CylindricalSolution,
    CylindricalSolutionCell,
    extract_cylindrical_solution,
)
from .formulas.boolean import is_false_expr, is_true_expr
from .normalization import normalize_formula, normalize_problem_variables


@dataclass(frozen=True)
class CADAdjacencyEdge:
    """Adjacency edge between two cylindrical CAD solution cells.

    The edge is semantic rather than purely syntactic: two selected cells are
    adjacent when their closures have a nonempty intersection that still lies
    in the original semialgebraic set. This avoids incorrectly connecting
    open components that merely touch in the ambient closure, such as
    ``x < 0`` and ``x > 0`` for the set ``x != 0``.
    """

    left: int
    right: int
    left_index: tuple[int, ...]
    right_index: tuple[int, ...]
    connector_formula: sp.Expr

    def as_pair(self) -> tuple[int, int]:
        return (self.left, self.right)


@dataclass(frozen=True)
class CADConnectedComponent:
    """A connected component of the selected CAD-cell adjacency graph."""

    index: int
    cells: tuple[CylindricalSolutionCell, ...]
    edges: tuple[CADAdjacencyEdge, ...]

    @property
    def dimension(self) -> int | None:
        if not self.cells:
            return None
        return max(cell.dimension for cell in self.cells)

    @property
    def bounded(self) -> bool | None:
        if not self.cells:
            return True
        return all(cell.bounded for cell in self.cells)

    def sample_point(self) -> dict[sp.Symbol, sp.Expr]:
        """Return a deterministic representative point for the component."""

        if not self.cells:
            return {}
        # Prefer a highest-dimensional cell so that samples usually lie in the
        # relative interior of a component rather than on its boundary.
        cell = max(self.cells, key=lambda item: (item.dimension, item.index))
        return cell.sample_point()

    @property
    def sample(self) -> Mapping[sp.Symbol, sp.Expr]:
        return self.sample_point()

    def as_formula(self, *, closed: bool = False) -> sp.Expr:
        if not self.cells:
            return sp.false
        pieces = [cell.as_formula(closed=closed) for cell in self.cells]
        return sp.Or(*pieces) if len(pieces) > 1 else pieces[0]


@dataclass(frozen=True)
class CADConnectivityGraph:
    """Roadmap-style adjacency graph for a CAD cell solution.

    This is not a full algebraic roadmap in the Canny/Basu-Pollack-Roy sense.
    It is a CAD-cell adjacency skeleton: nodes are selected cylindrical cells,
    edges record closure intersections that remain inside the solution set, and
    connected components of the graph give component-aware samples whenever the
    extracted CAD is sufficiently fine and faithful.
    """

    variables: tuple[sp.Symbol, ...]
    formula: sp.Expr
    cells: tuple[CylindricalSolutionCell, ...]
    edges: tuple[CADAdjacencyEdge, ...]
    components: tuple[CADConnectedComponent, ...]
    cylindrical_solution: CylindricalSolution | None = None

    @property
    def component_count(self) -> int:
        return len(self.components)

    @property
    def roadmap_nodes(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return tuple(cell.sample_point() for cell in self.cells)

    @property
    def roadmap_edges(self) -> tuple[tuple[int, int], ...]:
        return tuple(edge.as_pair() for edge in self.edges)

    @property
    def sample_points(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return tuple(component.sample_point() for component in self.components)


def _cell_formula(cell: CylindricalSolutionCell, *, closed: bool) -> sp.Expr:
    try:
        return cell.as_formula(closed=closed)
    except (ArithmeticError, NotImplementedError, sp.PolynomialError):
        parts = [level.as_formula(closed=closed) for level in cell.levels]
        return sp.And(*parts) if parts else sp.true


def _closures_intersect_inside_solution(
    left: CylindricalSolutionCell,
    right: CylindricalSolutionCell,
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
) -> tuple[bool, sp.Expr]:
    """Return whether two cell closures meet inside the solution formula."""

    connector = sp.And(_cell_formula(left, closed=True), _cell_formula(right, closed=True), formula)
    connector = sp.simplify(connector)
    if is_false_expr(connector):
        return False, sp.false
    try:
        from .decision import is_satisfiable

        return bool(is_satisfiable(connector, variables)), connector
    except (ArithmeticError, NotImplementedError, sp.PolynomialError):
        # Conservative fallback: only connect if SymPy can reduce the connector
        # to literal True. Do not guess adjacency from syntactic overlap.
        return bool(is_true_expr(connector)), connector


def _connected_components_from_edges(
    cells: Sequence[CylindricalSolutionCell],
    edges: Sequence[CADAdjacencyEdge],
) -> tuple[CADConnectedComponent, ...]:
    """Group CAD cells into connected components using union-find over certified adjacency edges."""
    n = len(cells)
    if n == 0:
        return ()
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    for edge in edges:
        union(edge.left, edge.right)

    groups: dict[int, list[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)

    out: list[CADConnectedComponent] = []
    for component_index, members in enumerate(sorted(groups.values(), key=lambda xs: min(xs))):
        member_set = set(members)
        component_edges = tuple(
            edge for edge in edges if edge.left in member_set and edge.right in member_set
        )
        out.append(
            CADConnectedComponent(
                index=component_index,
                cells=tuple(cells[i] for i in members),
                edges=component_edges,
            )
        )
    return tuple(out)


def _source_cad_edges(
    cylindrical_solution: CylindricalSolution,
    formula: sp.Expr,
) -> tuple[CADAdjacencyEdge, ...] | None:
    """Build adjacency from the certified source CAD when it is available.

    Structured CAD extraction retains the complete lifting tree. Its recursive
    closure-incidence test is much cheaper than launching a fresh satisfiability
    problem for every pair of already-certified cells, while preserving exact
    behavior at degenerating algebraic sections.
    """

    decomp = cylindrical_solution.source_decomposition
    if decomp is None or decomp.cad_result is None:
        return None

    result = decomp.cad_result
    cells = tuple(cell for cell in cylindrical_solution.cells if cell.selected)
    positions = {cell.index: pos for pos, cell in enumerate(cells)}
    source_cells = tuple(result.cell_set.cells)
    if set(positions) != {cell.index for cell in source_cells}:
        return None

    from .topology.components import component_cell_graph

    graph = component_cell_graph(source_cells, result.cad, cylindrical_solution.variables)
    edges: list[CADAdjacencyEdge] = []
    by_index = {cell.index: cell for cell in cells}
    for left_index in sorted(graph):
        left_pos = positions[left_index]
        for right_index in sorted(graph[left_index]):
            right_pos = positions[right_index]
            if right_pos <= left_pos:
                continue
            left = by_index[left_index]
            right = by_index[right_index]
            connector = sp.And(
                _cell_formula(left, closed=True),
                _cell_formula(right, closed=True),
                formula,
                evaluate=False,
            )
            edges.append(
                CADAdjacencyEdge(
                    left=left_pos,
                    right=right_pos,
                    left_index=left_index,
                    right_index=right_index,
                    connector_formula=connector,
                )
            )
    return tuple(edges)


def build_cad_adjacency_graph(
    cylindrical_solution: CylindricalSolution,
    *,
    formula: sp.Expr | None = None,
    max_pair_checks: int | None = None,
) -> CADConnectivityGraph:
    """Build a semantic adjacency graph on selected cylindrical CAD cells.

    ``max_pair_checks`` can be used by callers to avoid accidental quadratic
    explosions. When omitted, all pairs are checked. The implementation is
    exact for the pairs it checks because adjacency is certified by a real
    satisfiability query on the intersection of the two closed cell formulas
    and the original solution formula.
    """

    cells = tuple(cell for cell in cylindrical_solution.cells if cell.selected)
    variables = tuple(cylindrical_solution.variables)
    ambient_formula = normalize_formula(
        formula if formula is not None else cylindrical_solution.formula
    )
    source_edges = None
    if max_pair_checks is None:
        source_edges = _source_cad_edges(cylindrical_solution, ambient_formula)
    if source_edges is not None:
        components = _connected_components_from_edges(cells, source_edges)
        return CADConnectivityGraph(
            variables=variables,
            formula=ambient_formula,
            cells=cells,
            edges=source_edges,
            components=components,
            cylindrical_solution=cylindrical_solution,
        )

    edges: list[CADAdjacencyEdge] = []
    checks = 0
    for i, left in enumerate(cells):
        for j in range(i + 1, len(cells)):
            if max_pair_checks is not None and checks >= max_pair_checks:
                break
            checks += 1
            right = cells[j]
            adjacent, connector = _closures_intersect_inside_solution(
                left, right, ambient_formula, variables
            )
            if adjacent:
                edges.append(
                    CADAdjacencyEdge(
                        left=i,
                        right=j,
                        left_index=left.index,
                        right_index=right.index,
                        connector_formula=connector,
                    )
                )
        if max_pair_checks is not None and checks >= max_pair_checks:
            break
    components = _connected_components_from_edges(cells, edges)
    return CADConnectivityGraph(
        variables=variables,
        formula=ambient_formula,
        cells=cells,
        edges=tuple(edges),
        components=components,
        cylindrical_solution=cylindrical_solution,
    )


def _quadratic_zero_set_structurally_connected(
    factor: sp.Expr, variables: Sequence[sp.Symbol]
) -> bool | None:
    """Decide connectedness for affine or definite-quadratic zero sets cheaply."""
    vars_ = tuple(variables)
    try:
        poly = sp.Poly(factor, *vars_)
    except (sp.PolynomialError, TypeError, ValueError):
        return None
    degree = poly.total_degree()
    if degree <= 1:
        return True
    if degree != 2 or len(vars_) < 2:
        return None
    hessian = sp.ImmutableMatrix(sp.hessian(poly.as_expr(), vars_))
    if any(entry.free_symbols for entry in hessian):
        return None
    try:
        eigen_signs = []
        for value, mult in hessian.eigenvals().items():
            if value.is_positive is True:
                eigen_signs.extend([1] * int(mult))
            elif value.is_negative is True:
                eigen_signs.extend([-1] * int(mult))
            else:
                return None
        definite = len(eigen_signs) == len(vars_) and (
            all(s > 0 for s in eigen_signs) or all(s < 0 for s in eigen_signs)
        )
        if not definite:
            return None
        gradient = [sp.diff(poly.as_expr(), var) for var in vars_]
        critical = sp.solve(gradient, vars_, dict=True)
        if len(critical) != 1:
            return None
        extremum = sp.simplify(poly.as_expr().subs(critical[0]))
        if all(s > 0 for s in eigen_signs):
            if extremum.is_positive is True:
                return None  # empty zero set
            if extremum.is_nonpositive is True:
                return True
        else:
            if extremum.is_negative is True:
                return None  # empty zero set
            if extremum.is_nonnegative is True:
                return True
    except (TypeError, ValueError, NotImplementedError):
        return None
    return None


def factorized_equality_components(
    condition: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...] | None:
    """Return exact components for a simple reducible polynomial equality.

    The shortcut is used only when every squarefree factor zero set is itself
    connected. Inter-factor incidences are then decided from the original
    polynomial equations, avoiding reconstructed radical cell formulas. If any
    factor has more than one component, callers fall back to the general CAD
    topology path.
    """

    formula = normalize_formula(condition)
    if not isinstance(formula, sp.Equality):
        return None
    residual = sp.expand(formula.lhs - formula.rhs)
    try:
        poly = sp.Poly(residual, *variables)
        _, factors = sp.factor_list(poly.as_expr(), *variables)
    except (sp.PolynomialError, TypeError, ValueError, NotImplementedError):
        return None
    factor_exprs = tuple(sp.expand(factor) for factor, _mult in factors if factor.free_symbols)
    if len(factor_exprs) <= 1:
        return None

    factor_components: list[sp.Expr] = []
    for factor in factor_exprs:
        factor_formula = sp.Eq(factor, 0)
        structural = _quadratic_zero_set_structurally_connected(factor, variables)
        if structural is not True:
            graph = extract_cad_connectivity(factor_formula, variables)
            if graph.component_count != 1:
                return None
        factor_components.append(factor_formula)

    from .decision import is_satisfiable

    parent = list(range(len(factor_exprs)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        lroot, rroot = find(left), find(right)
        if lroot != rroot:
            parent[rroot] = lroot

    for left in range(len(factor_exprs)):
        for right in range(left + 1, len(factor_exprs)):
            intersection = sp.And(
                sp.Eq(factor_exprs[left], 0),
                sp.Eq(factor_exprs[right], 0),
            )
            if is_satisfiable(intersection, variables):
                union(left, right)

    groups: dict[int, list[sp.Expr]] = {}
    for index, factor_formula in enumerate(factor_components):
        groups.setdefault(find(index), []).append(factor_formula)
    return tuple(group[0] if len(group) == 1 else sp.Or(*group) for group in groups.values())


def extract_cad_connectivity(
    condition_or_solution: object,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    max_pair_checks: int | None = None,
) -> CADConnectivityGraph:
    """Extract a CAD adjacency/roadmap-style connectivity graph.

    The input may be either a formula-like condition or an existing
    ``CylindricalSolution``. Formula inputs are first decomposed by the current
    complete CAD backend and then converted into the cylindrical solution
    representation.
    """

    if isinstance(condition_or_solution, CylindricalSolution):
        cyl = condition_or_solution
        formula = cyl.formula
    else:
        formula = normalize_formula(condition_or_solution)
        vars_ = normalize_problem_variables(variables, formula)
        cyl = extract_cylindrical_solution(formula, vars_, selected_only=True)
    return build_cad_adjacency_graph(cyl, formula=formula, max_pair_checks=max_pair_checks)


__all__ = [
    "CADAdjacencyEdge",
    "CADConnectedComponent",
    "CADConnectivityGraph",
    "build_cad_adjacency_graph",
    "extract_cad_connectivity",
    "factorized_equality_components",
]
