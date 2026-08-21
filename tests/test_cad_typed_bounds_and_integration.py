from __future__ import annotations

import pytest
import sympy as sp

from semialg import (
    AlgebraicRootFunction,
    extract_cylindrical_solution,
    extract_explicit_cylindrical_solution,
    full_dimensional_cell_integral,
    intrinsic_cell_integral,
    verify_cad_cell_bounds,
)

pytestmark = pytest.mark.slow


def test_disk_uses_typed_algebraic_root_bounds_and_verifies():
    x, y = sp.symbols("x y", real=True)
    sol = extract_cylindrical_solution(sp.And(x >= 0, x <= 1, y**2 <= x), [x, y])
    cell = sol.full_dimensional_cells[0]
    y_level = cell.levels[1]
    assert isinstance(y_level.typed_lower, AlgebraicRootFunction)
    assert isinstance(y_level.typed_upper, AlgebraicRootFunction)
    assert sp.simplify(y_level.lower + sp.sqrt(x)) == 0
    assert sp.simplify(y_level.upper - sp.sqrt(x)) == 0
    assert verify_cad_cell_bounds(cell).verify()


def test_quintic_boundary_is_root_function_and_evaluates_exactly():
    x, y = sp.symbols("x y", real=True)
    sol = extract_cylindrical_solution(sp.And(x >= 0, x <= 1, y**5 + x * y + 1 <= 0), [x, y])
    cell = sol.full_dimensional_cells[0]
    upper = cell.levels[1].typed_upper
    assert isinstance(upper, AlgebraicRootFunction)
    assert upper.as_expr().func.__name__ == "root_of"
    value = upper.evaluate({x: 0})
    assert sp.simplify(value.as_expr() + 1) == 0
    assert verify_cad_cell_bounds(cell).verify()


def test_cubic_boundary_is_root_function_without_radical_requirement():
    x, y = sp.symbols("x y", real=True)
    sol = extract_cylindrical_solution(sp.And(x >= 0, x <= 1, y**3 + x * y + 1 <= 0), [x, y])
    cell = sol.full_dimensional_cells[0]
    upper = cell.levels[1].typed_upper
    assert isinstance(upper, AlgebraicRootFunction)
    assert sp.simplify(upper.evaluate({x: 0}).as_expr() + 1) == 0


def test_sphere_has_nested_root_bounds_in_y_and_z():
    x, y, z = sp.symbols("x y z", real=True)
    sol = extract_cylindrical_solution(x**2 + y**2 + z**2 <= 1, [x, y, z])
    cell = sol.full_dimensional_cells[0]
    assert isinstance(cell.levels[1].typed_lower, AlgebraicRootFunction)
    assert isinstance(cell.levels[2].typed_lower, AlgebraicRootFunction)
    assert cell.levels[2].typed_lower.base_variables == (x, y)
    assert sp.simplify(cell.levels[2].lower + sp.sqrt(1 - x**2 - y**2)) == 0
    assert verify_cad_cell_bounds(cell).verify()


def test_nested_cubic_z_bound_depends_on_x_and_y():
    x, y, z = sp.symbols("x y z", real=True)
    condition = sp.And(x >= 0, x <= 1, y >= 0, y <= 1, z**3 + x * z + y <= 0)
    sol = extract_cylindrical_solution(condition, [x, y, z])
    cell = sol.full_dimensional_cells[0]
    upper = cell.levels[2].typed_upper
    assert isinstance(upper, AlgebraicRootFunction)
    assert set(upper.base_variables) == {x, y}
    assert upper.polynomial.has(x, y, z)
    assert verify_cad_cell_bounds(cell).verify()


def test_unbounded_sector_preserves_infinity_and_openness():
    x, y = sp.symbols("x y", real=True)
    sol = extract_explicit_cylindrical_solution(sp.And(x > 0, y > x), [x, y])
    assert sol is not None
    cell = sol.full_dimensional_cells[0]
    assert cell.levels[0].upper == sp.oo
    assert cell.levels[0].lower_closed is False
    assert cell.levels[1].upper == sp.oo
    assert cell.levels[1].lower_closed is False
    assert verify_cad_cell_bounds(cell).verify()


def test_section_cell_preserves_equality_and_closedness():
    x, y = sp.symbols("x y", real=True)
    sol = extract_explicit_cylindrical_solution(sp.And(x >= 0, x <= 1, sp.Eq(y, x**2)), [x, y])
    assert sol is not None
    cell = sol.cells[0]
    section = cell.levels[1]
    assert section.is_section
    assert section.lower_closed and section.upper_closed
    assert sp.simplify(section.lower - x**2) == 0
    assert verify_cad_cell_bounds(cell).verify()


def test_mixed_dimensional_union_preserves_cell_dimensions():
    x, y = sp.symbols("x y", real=True)
    formula = sp.Or(
        sp.And(x > 0, x < 1, y > 0, y < 1),
        sp.And(sp.Eq(x, 2), sp.Eq(y, 0)),
    )
    sol = extract_cylindrical_solution(formula, [x, y])
    dimensions = {cell.dimension for cell in sol.cells}
    assert 2 in dimensions
    assert 0 in dimensions
    assert all(verify_cad_cell_bounds(cell).verify() for cell in sol.cells)


def test_root_order_swap_is_split_into_base_cells_with_certificates():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(x > -1, x < 1, y**2 <= x**2)
    sol = extract_cylindrical_solution(formula, [x, y])
    cells = sol.full_dimensional_cells
    assert len(cells) == 2
    left, right = cells
    assert left.levels[0].upper == 0
    assert right.levels[0].lower == 0
    assert left.levels[1].root_order is not None and left.levels[1].root_order.verify()
    assert right.levels[1].root_order is not None and right.levels[1].root_order.verify()
    assert sp.simplify(left.levels[1].lower - x) == 0
    assert sp.simplify(right.levels[1].lower + x) == 0


def test_root_function_evaluation_over_algebraic_base_cell():
    x, y = sp.symbols("x y", real=True)
    sol = extract_cylindrical_solution(sp.And(sp.Eq(x**2, 2), sp.Eq(y**3, x)), [x, y])
    positive = next(cell for cell in sol.cells if sp.N(cell.sample[x]) > 0)
    root = positive.levels[1].typed_lower
    assert isinstance(root, AlgebraicRootFunction)
    evaluated = root.evaluate({x: positive.sample[x]})
    assert sp.simplify(evaluated.as_expr() - 2 ** sp.Rational(1, 6)) == 0
    assert verify_cad_cell_bounds(positive).verify()


def test_full_dimensional_integration_adapter_handles_arbitrary_dimension():
    x, y, z = sp.symbols("x y z", real=True)
    sol = extract_explicit_cylindrical_solution(
        sp.And(x >= 0, x <= 1, y >= 0, y <= x, z >= 0, z <= x + y),
        [x, y, z],
    )
    assert sol is not None
    adapted = full_dimensional_cell_integral(sol.cells[0], 1, evaluate=True)
    assert adapted.certified_bounds
    assert sp.simplify(adapted.integral - sp.Rational(1, 2)) == 0


def test_intrinsic_integration_is_separate_and_uses_induced_metric():
    x, y, z = sp.symbols("x y z", real=True)
    sol = extract_explicit_cylindrical_solution(
        sp.And(x >= 0, x <= 1, y >= 0, y <= 1, sp.Eq(z, x + y)),
        [x, y, z],
    )
    assert sol is not None
    cell = sol.cells[0]
    assert cell.dimension == 2
    adapted = intrinsic_cell_integral(cell, 1, evaluate=True)
    assert adapted.intrinsic
    assert sp.simplify(adapted.metric_factor - sp.sqrt(3)) == 0
    assert sp.simplify(adapted.integral - sp.sqrt(3)) == 0
