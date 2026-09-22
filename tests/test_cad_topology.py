from __future__ import annotations

import pytest
import sympy as sp

from semialg.decomposition.components import components_from_cell_set
from semialg.decomposition.cylindrical import cad


def test_interval_boundary_interior_of_closure_exterior():
    x = sp.Symbol("x", real=True)
    interior = cad(x**2 <= 1, [x], operation="interior", return_result=True)
    closure = cad(x**2 < 1, [x], operation="closure", return_result=True)
    boundary = cad(x**2 <= 1, [x], operation="boundary", return_result=True)
    exterior = cad(x**2 <= 1, [x], operation="exterior", return_result=True)

    assert bool(interior.formula.subs({x: 0}))
    assert not bool(interior.formula.subs({x: 1}))
    assert bool(closure.formula.subs({x: -1}))
    assert bool(closure.formula.subs({x: 1}))
    assert bool(boundary.formula.subs({x: -1}))
    assert bool(boundary.formula.subs({x: 1}))
    assert not bool(boundary.formula.subs({x: 0}))
    assert bool(exterior.formula.subs({x: 2}))
    assert not bool(exterior.formula.subs({x: 1}))
    assert not bool(exterior.formula.subs({x: 0}))


@pytest.mark.slow
def test_disk_topological_operations():
    x, y = sp.symbols("x y", real=True)
    disk = x**2 + y**2 <= 1
    interior = cad(disk, [x, y], operation="interior", return_result=True)
    boundary = cad(disk, [x, y], operation="boundary", return_result=True)
    exterior = cad(disk, [x, y], operation="exterior", return_result=True)

    assert bool(interior.formula.subs({x: 0, y: 0}))
    assert not bool(interior.formula.subs({x: 1, y: 0}))
    assert bool(boundary.formula.subs({x: 1, y: 0}))
    assert bool(boundary.formula.subs({x: 0, y: 1}))
    assert not bool(boundary.formula.subs({x: 0, y: 0}))
    assert bool(exterior.formula.subs({x: 2, y: 0}))
    assert not bool(exterior.formula.subs({x: 1, y: 0}))


def test_components_from_cell_set_uses_closure_connected_cells():
    x = sp.Symbol("x", real=True)
    connected = cad(
        sp.Or((x >= -1) & (x <= 0), (x >= 0) & (x <= 1)), [x], output="cells", return_result=True
    )
    disconnected = cad(sp.Or(x < -1, x > 1), [x], output="cells", return_result=True)

    assert len(components_from_cell_set(connected.as_cell_set()).components) == 1
    assert len(components_from_cell_set(disconnected.as_cell_set()).components) == 2
