from __future__ import annotations

import sympy as sp

from semialg.decision.api import _prepare_solve_plan
from semialg.formula import parse_formula
from semialg.qe.complete import (
    _build_complete_qe_cad,
    _combine,
    _propagate_quantified_cell_truth,
    evaluate_formula_on_cell,
)
from semialg.solve.planner import affine_presolve, profile_semialgebraic_system


def test_planner_preparation_preserves_profile_and_affine_state():
    x, y = sp.symbols("x y", real=True)
    expr = sp.And(sp.Eq(x + y, 1), x >= 0, y >= 0)

    # Reference sequence copied from the pre-refactor planner: profile once,
    # then conditionally apply affine presolve and retain its working formula.
    profile = profile_semialgebraic_system(expr, (x, y))
    presolved = affine_presolve(expr, (x, y))
    expected_expr = presolved.formula if presolved.changed else expr
    expected_vars = presolved.variables if presolved.changed else (x, y)

    state = _prepare_solve_plan(expr, (x, y))
    assert state.profile == profile
    assert state.work_expr == expected_expr
    assert state.work_vars == expected_vars
    assert state.presolved == presolved
    assert state.reconstruct(sp.true) == sp.And(sp.Eq(x, 1 - y), sp.true)


def test_qe_truth_propagation_matches_direct_level_reduction():
    x, y = sp.symbols("x y", real=True)
    matrix = parse_formula(sp.And(x**2 <= 1, sp.Eq(y, x)))
    variables = (x, y)
    quantifiers = (("exists", y),)
    cad, _ = _build_complete_qe_cad(matrix, variables, quantifiers, allow_variety_cad=False)
    qmap = {y: "exists"}

    # Exact pre-refactor propagation loop retained here as a differential oracle.
    expected = {
        cell.index: evaluate_formula_on_cell(matrix, cell, variables)
        for cell in cad.cells_by_level[len(variables)]
    }
    for level in range(len(variables), 1, -1):
        var = variables[level - 1]
        grouped: dict[tuple[int, ...], list[bool]] = {}
        for idx, value in expected.items():
            grouped.setdefault(idx[:-1], []).append(value)
        expected = {parent: _combine(values, qmap[var]) for parent, values in grouped.items()}

    actual = _propagate_quantified_cell_truth(cad, matrix, variables, 1, qmap)
    assert actual == expected


def test_qe_backend_stage_preserves_collins_backend_choice():
    x, y = sp.symbols("x y", real=True)
    matrix = parse_formula(sp.And(x**2 <= 1, sp.Eq(y, x)))
    cad, notes = _build_complete_qe_cad(matrix, (x, y), (("exists", y),), allow_variety_cad=False)
    assert cad.backend != "groebner-variety"
    assert notes == ()
