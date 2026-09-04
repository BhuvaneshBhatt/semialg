from __future__ import annotations

import sympy as sp

from semialg import SemialgebraicRegion


def test_interval_complex_boundary_coboundary_duality_and_f_vector():
    x = sp.symbols("x", real=True)
    complex_ = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,)).cell_complex()

    assert complex_.f_vector == (2, 1)
    assert len(complex_.cells_by_dimension(0)) == 2
    assert len(complex_.cells_by_dimension(1)) == 1

    edge_pos = next(i for i, cell in enumerate(complex_.cells) if cell.dimension == 1)
    for vertex in complex_.boundary_indices(edge_pos):
        assert edge_pos in complex_.coboundary_indices(vertex)
        assert complex_.cells[vertex].dimension + 1 == complex_.cells[edge_pos].dimension


def test_cell_complex_euler_agrees_with_region_euler():
    x, y = sp.symbols("x y", real=True)
    region = SemialgebraicRegion(x**2 + y**2 <= 1, (x, y))
    complex_ = region.cell_complex()
    assert complex_.euler_characteristic() == region.euler_characteristic() == 1
