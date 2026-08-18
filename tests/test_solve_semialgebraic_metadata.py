import sympy as sp

from semialg import solve_semialgebraic
from semialg.decision import SemialgebraicSolution


def test_solve_semialgebraic_returns_metadata_for_interval():
    x = sp.symbols("x", real=True)
    sol = solve_semialgebraic([x >= 0, x <= 1], [x])
    assert isinstance(sol, SemialgebraicSolution)
    assert sol.satisfiable is True
    assert sol.nonempty is True
    assert sol.empty is False
    assert sol.variables == (x,)
    assert sol.dimension == 1
    assert sol.bounded is True
    assert sol.closed is True
    assert sol.compact is True
    assert sol.simplified_constraints
    assert sol.diagnostics["solver_stage"] == "semialgebraic_solution"


def test_solve_semialgebraic_metadata_for_unsat_system():
    x = sp.symbols("x", real=True)
    sol = solve_semialgebraic([x > 0, x <= 0], [x])
    assert sol.satisfiable is False
    assert sol.empty is True
    assert sol.formula == sp.false
    assert sol.samples == ()


def test_solve_semialgebraic_parameter_conditions_are_recorded():
    x, a, b = sp.symbols("x a b", real=True)
    sol = solve_semialgebraic([sp.Eq(x**2 + a * x + b, 0)], [x], parameters=[a, b], count=0)
    assert sol.variables == (x,)
    assert sol.parameters == (a, b)
    assert sol.parameter_conditions == (a**2 - 4 * b >= 0)
    assert sol.diagnostics["has_parameter_conditions"] is True


def test_solve_semialgebraic_records_vertical_cells_for_supported_2d_slice():
    x, y = sp.symbols("x y", real=True)
    sol = solve_semialgebraic([x >= 0, x <= 1, y >= x, y <= 1], [x, y], count=0)
    assert sol.satisfiable is True
    assert sol.dimension == 2
    assert sol.bounded is True
    assert len(sol.cells) == 1
    cell = sol.cells[0]
    assert cell.x_interval == (sp.Integer(0), sp.Integer(1))
    assert cell.y_bounds == ((x, sp.Integer(1)),)
