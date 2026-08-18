from __future__ import annotations

import sympy as sp

from semialg import component_instances
from semialg.decomposition.components import cell_adjacency_graph, components_from_cell_set
from semialg.decomposition.cylindrical import cad


def test_component_inst_01():
    x = sp.Symbol("x", real=True)
    result = component_instances(sp.Or(x < -1, x > 1), [x])
    assert result.status == "complete"
    assert len(result.components) == 2
    assert all(x in sample for sample in result.instances)


def test_component_inst_02():
    x = sp.Symbol("x", real=True)
    cad_result = cad(sp.Or(x < -1, x > 1), [x], output="cells")
    graph = cell_adjacency_graph(cad_result.as_cell_set())
    assert len(graph.nodes) == len(cad_result.cells)
    assert not any(graph.neighbors(node) for node in graph.nodes)


def test_component_inst_03():
    x = sp.Symbol("x", real=True)
    cad_result = cad(x**2 <= 1, [x], output="cells")
    result = components_from_cell_set(cad_result.as_cell_set())
    assert len(result.components) == 1
    assert result.components[0].dimension == 1
