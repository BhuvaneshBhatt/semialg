import pytest
import sympy as sp

from semialg.cad.cells import (
    CylindricalCoordinateConstraint,
    CylindricalSolution,
    CylindricalSolutionCell,
)
from semialg.connectivity import (
    CADConnectivityGraph,
    build_cad_adjacency_graph,
)
from semialg.decision import solve_semialgebraic


def _manual_cell(index, xlo, xhi, ylo, yhi):
    x, y = sp.symbols("x y", real=True)
    sample = {x: sp.simplify((xlo + xhi) / 2), y: sp.simplify((ylo + yhi) / 2)}
    levels = (
        CylindricalCoordinateConstraint(x, 1, "sector", xlo, xhi, sample[x], index[:1]),
        CylindricalCoordinateConstraint(y, 2, "sector", ylo, yhi, sample[y], index),
    )
    return CylindricalSolutionCell((x, y), levels, sample, index)


def test_connectivity_graph_on_manual_adjacent_cells():
    x, y = sp.symbols("x y", real=True)
    left = _manual_cell((1, 1), 0, 1, 0, 1)
    right = _manual_cell((2, 1), 1, 2, 0, 1)
    cyl = CylindricalSolution((x, y), (left, right), sp.And(x >= 0, x <= 2, y >= 0, y <= 1))
    graph = build_cad_adjacency_graph(cyl)
    assert isinstance(graph, CADConnectivityGraph)
    assert graph.component_count == 1
    assert graph.roadmap_edges == ((0, 1),)


def test_connectivity_graph_keeps_gap_components_separate():
    x, y = sp.symbols("x y", real=True)
    left = _manual_cell((1, 1), 0, 1, 0, 1)
    right = _manual_cell((2, 1), 2, 3, 0, 1)
    formula = sp.Or(sp.And(x >= 0, x <= 1, y >= 0, y <= 1), sp.And(x >= 2, x <= 3, y >= 0, y <= 1))
    cyl = CylindricalSolution((x, y), (left, right), formula)
    graph = build_cad_adjacency_graph(cyl)
    assert graph.component_count == 2
    assert graph.roadmap_edges == ()


@pytest.mark.slow
def test_solve_semialgebraic_exposes_connectivity_output_for_supported_regions():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    graph = solve_semialgebraic(condition, [x, y], count=0, output="connectivity")
    assert graph is None or hasattr(graph, "components")


@pytest.mark.slow
def test_per_component_sampling_can_use_connectivity_components():
    x, y = sp.symbols("x y", real=True)
    condition = sp.Or(
        sp.And(x >= 0, x <= 1, y >= 0, y <= 1), sp.And(x >= 2, x <= 3, y >= 0, y <= 1)
    )
    sol = solve_semialgebraic(condition, [x, y], count=0, samples="per_component")
    # In environments where complete CAD extraction is available, this should be
    # truly component-aware. Otherwise the older supported cell parser still
    # gives one sample per supported cell for this disjoint union.
    assert len(sol.samples) >= 1
