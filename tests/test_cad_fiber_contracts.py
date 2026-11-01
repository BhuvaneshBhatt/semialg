"""Exact projection, reconstruction and root-certificate contracts."""

from dataclasses import replace
from itertools import product

import pytest
import sympy as sp

from semialg import cad
from semialg.cad_algorithms.bounds import verify_cad_cell_bounds
from semialg.cad_algorithms.cells import extract_cylindrical_solution
from semialg.reconstruct.merge import compressed_formula_from_cells, dnf_formula_from_cells
from semialg.topology.operations import cells_formula


@pytest.mark.parametrize("method", ["project", "project_cell_set"])
@pytest.mark.parametrize("family", ["diagonal", "grid", "algebraic", "empty"])
def test_finite_projection(method, family):
    x, y = sp.symbols("x y", real=True)
    if family == "diagonal":
        expression = sp.And(sp.Eq(y * y, 1), sp.Eq(x, y))
    elif family == "grid":
        expression = sp.And(sp.Eq(x * x, 1), sp.Eq(y * y, 1))
    elif family == "algebraic":
        expression = sp.And(sp.Eq(x * x, 2), sp.Eq(x, y))
    else:
        expression = sp.And(sp.Eq(x * x + 1, 0, evaluate=False), sp.Eq(y, 0))
    source = cad(expression, (x, y), output="function")
    assert source.cad.tower.metadata["variety_only"]
    result = getattr(source, method)((x,))
    for point in (-2, -sp.sqrt(2), -1, 0, 1, sp.sqrt(2), 2):
        expected = (
            bool(sp.Eq(point**2, 2))
            if family == "algebraic"
            else (point in (-1, 1) and family != "empty")
        )
        actual = result(point) if method == "project" else bool(result.formula.subs(x, point))
        assert actual is expected


@pytest.mark.parametrize("method", ["project", "project_cell_set"])
def test_projection_preserves_point_correlations(method):
    x, y, z = sp.symbols("x y z", real=True)
    expression = sp.And(sp.Eq(y * y, 1), sp.Eq(x, y), sp.Eq(z, 0))
    source = cad(expression, (x, y, z), output="function")
    result = getattr(source, method)((x, y))
    for first, second in product((-1, 0, 1), repeat=2):
        expected = first == second and first in (-1, 1)
        actual = (
            result(first, second)
            if method == "project"
            else bool(result.formula.subs({x: first, y: second}))
        )
        assert actual is expected


@pytest.mark.parametrize("method", ["project", "project_cell_set"])
def test_ambient_projection_keeps_intervals(method):
    x, y = sp.symbols("x y", real=True)
    source = cad(sp.Eq(x * x + y * y, 1), (x, y), output="function")
    result = getattr(source, method)((x,))
    for point in (-2, -1, -sp.Rational(1, 2), 0, sp.Rational(1, 2), 1, 2):
        expected = -1 <= point <= 1
        actual = result(point) if method == "project" else bool(result.formula.subs(x, point))
        assert actual is bool(expected)


@pytest.mark.parametrize("left_closed,right_closed", list(product((False, True), repeat=2)))
@pytest.mark.parametrize("dimension", [1, 2])
def test_closed_compression(left_closed, right_closed, dimension):
    x, y = sp.symbols("x y", real=True)
    variables = (x, y)[:dimension]
    interval = sp.And(x >= 0 if left_closed else x > 0, x <= 1 if right_closed else x < 1)
    expression = interval if dimension == 1 else sp.And(interval, y >= 0, y <= 1)
    source = cad(expression, variables, strategy="full", return_result=True)
    compressed = compressed_formula_from_cells(
        source.cells, variables, source.cad.cells_by_level, closed=True
    )
    expanded = dnf_formula_from_cells(
        source.cells, variables, source.cad.cells_by_level, closed=True
    )
    for values in product((-1, 0, sp.Rational(1, 2), 1, 2), repeat=dimension):
        point = dict(zip(variables, values, strict=True))
        expected = all(0 <= value <= 1 for value in values)
        assert bool(compressed.formula.subs(point)) is expected
        assert bool(expanded.subs(point)) is expected


def test_topology_reconstruction_closes_merged_bounds():
    x = sp.Symbol("x", real=True)
    source = cad(sp.And(x > 0, x <= 1), (x,), return_result=True)
    result = cells_formula(source.cells, (x,), source.cad.cells_by_level, closed=True)
    assert bool(result.subs(x, 0))
    assert bool(result.subs(x, 1))
    assert not bool(result.subs(x, -1))


@pytest.fixture
def root_cell():
    x, y = sp.symbols("x y", real=True)
    source = extract_cylindrical_solution(sp.And(x > 0, x * y * y <= 1), (x, y))
    return source.full_dimensional_cells[0]


@pytest.mark.parametrize(
    "field",
    [
        "polynomial",
        "fiber_variable",
        "root_index",
        "base_variables",
        "base_index",
        "stack_root_index",
    ],
)
def test_root_certificate_binding(root_cell, field):
    x, y = root_cell.variables
    level = root_cell.levels[1]
    lower = level.typed_lower
    changes = {
        "polynomial": lower.polynomial + y,
        "fiber_variable": x,
        "root_index": 1,
        "base_variables": (),
        "base_index": (0,),
        "stack_root_index": 1,
    }
    wrong = replace(lower, certificate=replace(lower.certificate, **{field: changes[field]}))
    changed = replace(root_cell, levels=(root_cell.levels[0], replace(level, lower_bound=wrong)))
    assert not verify_cad_cell_bounds(changed).verify()
    assert not wrong.is_regular
    assert wrong.as_expr().func.__name__ == "root_of"
    assert wrong.as_expr().subs(x, 1) == -1
    assert not wrong.compare_certified(level.typed_upper).certified
    with pytest.raises(ValueError):
        wrong.derivative_expr(x)


def test_root_copy_does_not_reuse_another_branch(root_cell):
    x, _ = root_cell.variables
    lower = root_cell.levels[1].typed_lower
    changed = replace(lower, root_index=1)
    assert changed.as_expr().subs(x, 1) == changed.evaluate({x: 1}).as_expr() == 1
    assert not changed.is_regular


@pytest.mark.parametrize(
    "field", ["sign_invariant", "stack_order_verified", "sample_root_verified"]
)
def test_unproved_radical_branch_is_opaque(root_cell, field):
    x, _ = root_cell.variables
    lower = root_cell.levels[1].typed_lower
    changed = replace(lower, certificate=replace(lower.certificate, **{field: False}))
    assert changed.as_expr().func.__name__ == "root_of"
    assert changed.as_expr().subs(x, 1) == -1


@pytest.mark.parametrize("value", [None, -1, 0, 1])
def test_root_comparison_uses_requested_point(root_cell, value):
    x, _ = root_cell.variables
    level = root_cell.levels[1]
    point = {} if value is None else {x: value}
    result = level.typed_lower.compare_certified(level.typed_upper, base_point=point)
    if value == 1:
        assert result.certified and result.relation == -1
        assert result.scope == "point"
    else:
        assert not result.certified
        assert result.relation is None


def test_matching_root_certificate_remains_usable(root_cell):
    level = root_cell.levels[1]
    lower, upper = level.typed_lower, level.typed_upper
    result = lower.compare_certified(upper)
    assert result.certified and result.relation == -1
    assert result.scope == "base-cell"
    assert lower.is_regular and upper.is_regular
    assert verify_cad_cell_bounds(root_cell).verify()


@pytest.mark.parametrize("relation", [sp.Eq, sp.Ne, sp.Lt, sp.Le, sp.Gt, sp.Ge])
@pytest.mark.parametrize("interface", ["eliminate", "planner"])
def test_virtual_substitution_handles_zero_atoms(relation, interface):
    from semialg.qe.virtual_substitution import (
        can_use_quadratic_vs,
        eliminate_quadratic_variable,
        try_quadratic_vs_qe,
    )

    x = sp.Symbol("x", real=True)
    matrix = sp.And(relation(0, 0, evaluate=False), x * x < 1, evaluate=False)
    assert can_use_quadratic_vs(matrix, x)
    if interface == "eliminate":
        result = eliminate_quadratic_variable(matrix, x)
    else:
        result = try_quadratic_vs_qe((x,), (("exists", x),), matrix).formula
    assert result is relation(0, 0)
