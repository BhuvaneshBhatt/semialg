from __future__ import annotations

import sympy as sp

from semialg import (
    decompose_cylindrical_formula_to_vertical_bounds_2d,
    decompose_implicit_formula,
    extract_symbolic_box_bounds,
    integrate_over_region,
    reduce_region_integral,
    region_boundary,
    semialgebraic_level_function,
)


def test_semialgebraic_level_function_for_boolean_inequalities():
    x, y = sp.symbols("x y", real=True)
    level = semialgebraic_level_function(sp.And(x**2 + y**2 <= 1, y >= 0), [x, y])
    assert sp.simplify(level - sp.Max(x**2 + y**2 - 1, -y)) == 0


def test_decompose_implicit_formula_splits_dnf_and_ignores_unequalities():
    x, y = sp.symbols("x y", real=True)
    pieces = decompose_implicit_formula(sp.Or(x < 0, sp.And(x >= 1, sp.Eq(y, 0), x != 2)), [x, y])
    formulas = {sp.sstr(piece.as_formula()) for piece in pieces}
    assert sp.sstr(x <= 0) in formulas
    assert sp.sstr(sp.And(1 - x <= 0, sp.Eq(y, 0))) in formulas


def test_extract_symbolic_box_bounds():
    x, y = sp.symbols("x y", real=True)
    box = extract_symbolic_box_bounds(sp.And(x >= 0, x <= 1, y >= -2, y <= 3), [x, y])
    assert box.limits == ((x, 0, 1), (y, -2, 3))


def test_cylindrical_vertical_bounds_reduce_triangle_like_cell():
    x, y = sp.symbols("x y", real=True)
    cells = decompose_cylindrical_formula_to_vertical_bounds_2d(
        sp.And(x >= 0, x <= 1, y >= x, y <= 1), [x, y]
    )
    full = [cell for cell in cells if cell.is_full_dimensional]
    assert len(full) == 1
    assert full[0].x_interval == (0, 1)
    assert full[0].y_bounds == ((x, 1),)
    assert integrate_over_region(1, sp.And(x >= 0, x <= 1, y >= x, y <= 1), [x, y]) == sp.Rational(
        1, 2
    )


def test_cylindrical_vertical_bounds_handle_disjoint_union_in_integral_reduction():
    x, y = sp.symbols("x y", real=True)
    cond = sp.Or(
        sp.And(x >= 0, x <= 1, y >= 0, y <= 1),
        sp.And(x >= 2, x <= 3, y >= 0, y <= 1),
    )
    reduced = reduce_region_integral(x + y, cond, [x, y])
    assert len(reduced.pieces) == 2
    assert integrate_over_region(x + y, cond, [x, y]) == 4


def test_region_boundary_uses_vertical_bounds_for_cylindrical_region():
    x, y = sp.symbols("x y", real=True)
    boundary = region_boundary(sp.And(x >= 0, x <= 1, y >= x, y <= 1), [x, y])
    assert boundary.subs({x: sp.Rational(1, 2), y: sp.Rational(1, 2)}) == sp.true
    assert boundary.subs({x: sp.Rational(1, 2), y: 1}) == sp.true
    assert boundary.subs({x: sp.Rational(1, 2), y: sp.Rational(3, 4)}) == sp.false
