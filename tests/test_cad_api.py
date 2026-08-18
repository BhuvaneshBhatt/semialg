from __future__ import annotations

import pytest
import sympy as sp

from semialg.decomposition import CADFunction, CellSet, cad, cad_text

pytestmark = pytest.mark.slow


def test_cad_api_01():
    x = sp.Symbol("x", real=True)
    result = cad(x**2 - 1 >= 0, [x])
    assert result.status == "complete"
    assert len(result.cells) >= 2
    assert bool(result.formula.subs({x: -2}))
    assert not bool(result.formula.subs({x: 0}))


def test_cad_api_02():
    cell_set = cad_text("x^2 <= 1", variables=["x"], output="cells", return_result=False)
    assert isinstance(cell_set, CellSet)
    assert len(cell_set) > 0
    assert cell_set.sample_points()


def test_cad_api_03():
    x, y = sp.symbols("x y", real=True)
    fn = cad((x >= 0) & (y >= 0), [x, y], output="function", return_result=False)
    assert isinstance(fn, CADFunction)
    assert fn.contains({x: 1, y: 1})
    assert not fn.contains({x: -1, y: 1})


def test_cad_api_04():
    x = sp.Symbol("x", real=True)
    interior = cad(x**2 <= 1, [x], operation="interior")
    closure = cad(x**2 < 1, [x], operation="closure")
    assert bool(interior.formula.subs({x: 0}))
    assert not bool(interior.formula.subs({x: 1}))
    assert bool(closure.formula.subs({x: 1}))


def test_cad_api_05():
    import semialg

    assert semialg.cad is not None


def test_cad_api_06():
    x, y = sp.symbols("x y", real=True)
    fn = cad((x >= 0) & (y >= 0), [x, y], output="function", return_result=False)
    projected = fn.project([x])
    assert isinstance(projected, CADFunction)
    assert projected.contains({x: 1})
    projected_cells = fn.project_cell_set([x])
    assert isinstance(projected_cells, CellSet)
    assert bool(projected_cells.to_formula().subs({x: 1}))
    restricted = fn.restrict({x: 1})
    assert isinstance(restricted, CADFunction)
    assert restricted.contains({y: 1})
    assert not restricted.contains({y: -1})
    assert fn.tree_by_index()
    tree = cad((x >= 0) & (y >= 0), [x, y], output="tree", return_result=False)
    assert tree
