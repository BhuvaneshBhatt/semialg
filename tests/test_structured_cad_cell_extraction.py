import pytest
import sympy as sp

from semialg import (
    extract_structured_cad_cells,
    extract_vertical_bounds_from_cad_2d,
    integrate_over_region,
    reduce_region_integral,
    region_boundary,
    solve_semialgebraic,
)

pytestmark = pytest.mark.slow


def test_extract_structured_cad_cells_for_parabolic_region():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x >= 0, x <= 1, y**2 <= x)
    decomposition = extract_structured_cad_cells(condition, [x, y])
    assert decomposition.variables == (x, y)
    assert decomposition.cells
    assert any(cell.dimension == 2 for cell in decomposition.cells)
    assert all(set(cell.sample).issubset({x, y}) for cell in decomposition.cells)


def test_extract_vertical_bounds_from_complete_cad_for_nonlinear_stack():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x >= 0, x <= 1, y**2 <= x)
    cells = extract_vertical_bounds_from_cad_2d(condition, [x, y], full_dimensional_only=True)
    assert cells
    assert any(cell.x_interval[0] == 0 and cell.x_interval[1] == 1 for cell in cells)
    formulas = [cell.as_formula(use_source=False) for cell in cells]
    assert any(formula.has(sp.sqrt(x)) for formula in formulas)


def test_reduce_region_integral_uses_complete_cad_vertical_bounds():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x >= 0, x <= 1, y**2 <= x)
    reduced = reduce_region_integral(1, condition, [x, y])
    assert reduced.method == "complete_cad_vertical_bounds_2d"
    assert reduced.pieces
    assert sp.simplify(reduced.unevaluated_sum().doit() - sp.Rational(4, 3)) == 0


def test_integrate_over_region_complete_cad_vertical_bounds():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x >= 0, x <= 1, y**2 <= x)
    assert sp.simplify(integrate_over_region(1, condition, [x, y]) - sp.Rational(4, 3)) == 0


def test_region_boundary_uses_complete_cad_bounds_for_parabolic_region():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x >= 0, x <= 1, y**2 <= x)
    boundary = region_boundary(condition, [x, y])
    assert boundary.has(sp.sqrt(x))
    assert bool(boundary.subs({x: sp.Rational(1, 4), y: sp.Rational(1, 2)}))


def test_solve_semialgebraic_cells_uses_complete_cad_fallback():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x >= 0, x <= 1, y**2 <= x)
    solution = solve_semialgebraic(condition, [x, y], count=0, output="cells")
    assert solution.cells
    assert any(cell.sample_point() for cell in solution.cells)
