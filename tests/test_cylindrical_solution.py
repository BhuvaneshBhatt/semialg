import pytest
import sympy as sp

from semialg import (
    CylindricalSolution,
    extract_cylindrical_solution,
    solve_semialgebraic,
)

pytestmark = pytest.mark.slow


def test_extract_cylindrical_solution_for_vertical_region_has_nested_bounds():
    x, y = sp.symbols("x y", real=True)
    sol = extract_cylindrical_solution(sp.And(x >= 0, x <= 1, y >= x, y <= 1), [x, y])
    assert isinstance(sol, CylindricalSolution)
    assert sol.variables == (x, y)
    assert sol.cells
    full = sol.full_dimensional_cells
    assert full
    cell = full[0]
    assert len(cell.levels) == 2
    assert cell.levels[0].variable == x
    assert cell.levels[1].variable == y
    assert cell.levels[1].lower.has(x) or cell.levels[1].upper.has(x)
    assert cell.sample_point()[x] >= 0


def test_solve_semialgebraic_output_cylindrical_returns_solution():
    x, y = sp.symbols("x y", real=True)
    sol = solve_semialgebraic([x >= 0, x <= 1, y >= x, y <= 1], [x, y], count=0)
    assert sol.cylindrical_solution is not None
    cyl = solve_semialgebraic(
        [x >= 0, x <= 1, y >= x, y <= 1], [x, y], count=0, output="cylindrical"
    )
    assert isinstance(cyl, CylindricalSolution)
    assert cyl.full_dimensional_cells


def test_cylindrical_cell_formula_and_limits():
    x, y = sp.symbols("x y", real=True)
    cyl = extract_cylindrical_solution(sp.And(x >= 0, x <= 1, y >= x, y <= 1), [x, y])
    cell = cyl.full_dimensional_cells[0]
    assert cell.as_formula() != sp.false
    limits = cell.iterated_limits()
    assert limits[0][0] == x
    assert limits[1][0] == y
