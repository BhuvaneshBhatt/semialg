from __future__ import annotations

import sympy as sp

from semialg import SemialgebraicRegion
from semialg.cad_algorithms.cells import extract_structured_cad_cells
from semialg.cad_region import triangulate_cad_cell, triangulate_cad_cells


def _adjacent_rectangle_cells():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(
        x >= 0,
        x <= 2,
        y >= 0,
        y <= 1,
        sp.Or(x <= 1, x >= 1, evaluate=False),
        evaluate=False,
    )
    cad = SemialgebraicRegion(formula, (x, y)).as_cad_region()
    decomp = extract_structured_cad_cells(cad.result)
    return tuple(cell for cell in decomp.full_dimensional_cells if cell.bounded)


def test_single_cell_mesh_has_valid_simplices_and_vertices():
    cells = _adjacent_rectangle_cells()
    mesh = triangulate_cad_cell(cells[0], slices=2, inset=0)
    assert mesh.dimension == 2
    assert all(len(simplex) == 3 for simplex in mesh.simplices)
    assert all(len(set(simplex)) == 3 for simplex in mesh.simplices)
    assert all(0 <= index < len(mesh.coordinates) for s in mesh.simplices for index in s)


def test_multicell_mesh_shares_boundary_vertices_and_is_conforming():
    cells = _adjacent_rectangle_cells()
    assert len(cells) >= 2
    mesh = triangulate_cad_cells(cells, slices=3, require_conforming=True)
    assert mesh.conforming is True
    assert mesh.shared_vertex_count > 0
    assert len(mesh.simplex_cells) == len(mesh.simplices)
    assert any(len(source_cells) > 1 for source_cells in mesh.vertex_cells)


def test_mesh_coordinates_stay_inside_rectangle_closure():
    cells = _adjacent_rectangle_cells()
    mesh = triangulate_cad_cells(cells, slices=2, require_conforming=True)
    for x_value, y_value in mesh.coordinates:
        assert -1e-12 <= float(x_value) <= 2 + 1e-12
        assert -1e-12 <= float(y_value) <= 1 + 1e-12
