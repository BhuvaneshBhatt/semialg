import pytest
import sympy as sp

from semialg import VerticalBoundCell2D, solve_semialgebraic

pytestmark = pytest.mark.slow


def test_output_cells_for_single_vertical_slice_region():
    x, y = sp.symbols("x y", real=True)

    cells = solve_semialgebraic(
        [x >= 0, x <= 1, y >= x, y <= 1],
        [x, y],
        count=0,
        output="cells",
    )

    assert len(cells) == 1
    cell = cells[0]
    assert isinstance(cell, VerticalBoundCell2D)
    assert cell.dimension == 2
    assert cell.bounded is True
    assert cell.is_full_dimensional is True
    assert cell.x_interval == (sp.Integer(0), sp.Integer(1))
    assert cell.y_bounds == ((x, sp.Integer(1)),)
    assert cell.sample_point() == {x: sp.Rational(1, 2), y: sp.Rational(3, 4)}


def test_result_records_cells_for_disjoint_vertical_slice_union():
    x, y = sp.symbols("x y", real=True)
    left_square = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    right_square = sp.And(x >= 2, x <= 3, y >= 0, y <= 1)

    sol = solve_semialgebraic(sp.Or(left_square, right_square), [x, y], count=0)

    assert sol.satisfiable is True
    assert sol.dimension == 2
    assert sol.bounded is True
    assert sol.closed is True
    assert sol.compact is True
    assert len(sol.cells) == 2
    ordered = sorted(sol.cells, key=lambda cell: sp.default_sort_key(cell.x_interval[0]))
    assert [cell.x_interval for cell in ordered] == [
        (sp.Integer(0), sp.Integer(1)),
        (sp.Integer(2), sp.Integer(3)),
    ]
    assert [cell.sample_point() for cell in ordered] == [
        {x: sp.Rational(1, 2), y: sp.Rational(1, 2)},
        {x: sp.Rational(5, 2), y: sp.Rational(1, 2)},
    ]


def test_cells_metadata_flag_is_present():
    x, y = sp.symbols("x y", real=True)

    sol = solve_semialgebraic([x >= 0, x <= 1, y >= 0, y <= 1], [x, y], count=0)

    assert sol.diagnostics["solver_stage"] == "semialgebraic_solution"
    assert sol.diagnostics["capabilities"]["vertical_cells"] == "two_dimensional"


def test_open_vertical_slice_has_cells_but_is_not_closed_or_compact():
    x, y = sp.symbols("x y", real=True)

    sol = solve_semialgebraic([x > 0, x < 1, y > 0, y < 1], [x, y], count=0)

    assert len(sol.cells) == 1
    assert sol.bounded is True
    assert sol.closed is False
    assert sol.compact is False
