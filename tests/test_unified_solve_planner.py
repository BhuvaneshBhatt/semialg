from __future__ import annotations

import pytest
import sympy as sp

from semialg import is_satisfiable, solve_semialgebraic


def _methods(result):
    return tuple(step["method"] for step in result.diagnostics.get("planner_steps", ()))


def test_affine_presolve_reconstructs_eliminated_output_variable():
    x, y = sp.symbols("x y", real=True)
    result = solve_semialgebraic(sp.And(sp.Eq(x + y, 1), x >= 0, y >= 0), (x, y), count=0)

    assert result.satisfiable
    assert "affine_presolve" in _methods(result)
    assert sp.Eq(x, 1 - y) in sp.And.make_args(result.formula)
    assert result.contains({x: sp.Rational(1, 4), y: sp.Rational(3, 4)})
    assert not result.contains({x: 1, y: 1})


def test_coupled_linear_system_uses_fourier_motzkin_not_cad():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(x >= 0, y >= 0, x + y <= 1)
    result = solve_semialgebraic(formula, (x, y), count=0)

    assert result.method == "linear_fourier_motzkin"
    assert result.satisfiable
    assert _methods(result)[-1] == "linear_fourier_motzkin"
    assert "cad_qe" not in _methods(result)


def test_linear_infeasibility_is_proved_without_cad():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(x >= 0, y >= 0, x + y <= -1)
    result = solve_semialgebraic(formula, (x, y), count=0)

    assert result.method == "linear_fourier_motzkin"
    assert result.satisfiable is False
    assert result.formula is sp.false


def test_independent_variable_blocks_are_solved_separately():
    x, y = sp.symbols("x y", real=True)
    result = solve_semialgebraic(sp.And(x**2 <= 1, y**2 >= 4), (x, y), count=0)

    assert result.satisfiable
    assert result.method.startswith("incidence_decomposition[")
    assert _methods(result)[-1] == "incidence_decomposition"
    assert len(result.diagnostics["child_plans"]) == 2
    assert bool(result.formula.subs({x: 0, y: 2})) is True
    assert bool(result.formula.subs({x: 2, y: 2})) is False


def test_small_boolean_union_is_solved_branchwise():
    x, y = sp.symbols("x y", real=True)
    formula = sp.Or(
        sp.And(x >= 0, y >= 0, x + y <= 1),
        sp.And(x >= 2, y >= 0, x + y <= 3),
    )
    result = solve_semialgebraic(formula, (x, y), count=0)

    assert result.satisfiable
    assert _methods(result)[-1] == "boolean_branch_decomposition"
    assert len(result.diagnostics["branch_plans"]) == 2
    assert bool(result.formula.subs({x: 0, y: 0})) is True
    assert bool(result.formula.subs({x: 2, y: 0})) is True
    assert bool(result.formula.subs({x: sp.Rational(3, 2), y: 0})) is False


def test_affine_presolve_keeps_rur_provenance_for_finite_system():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x**2 + y**2, 1), sp.Eq(x, y))
    result = solve_semialgebraic(formula, (x, y), count=10)

    assert result.method == "rational_univariate"
    assert _methods(result)[0] == "affine_presolve"
    assert _methods(result)[-1] == "rational_univariate"
    assert len(result.samples) == 2


def test_positive_dimensional_equality_uses_groebner_presolve_before_cad():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x * y, 0), x >= 0)
    result = solve_semialgebraic(formula, (x, y), count=0)

    methods = _methods(result)
    assert "groebner_presolve" in methods
    assert methods[-1] == "cad_qe"
    assert bool(result.formula.subs({x: 0, y: 5})) is True
    assert bool(result.formula.subs({x: 2, y: 0})) is True
    assert bool(result.formula.subs({x: -1, y: 0})) is False


def test_explicit_linear_method_rejects_nonlinear_system():
    x = sp.symbols("x", real=True)
    with pytest.raises(NotImplementedError):
        solve_semialgebraic(x**2 <= 1, (x,), method="linear", count=0)


def test_is_satisfiable_reuses_linear_fourier_motzkin_path():
    x, y, z = sp.symbols("x y z", real=True)
    formula = sp.And(x >= 0, y >= 0, z >= 0, x + y + z <= -1)
    result = is_satisfiable(formula, (x, y, z), return_result=True)

    assert result.satisfiable is False
    assert result.method == "linear_fourier_motzkin"


def test_exact_ideal_dimension_rejects_dependent_equations_as_finite():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x * y, 0), sp.Eq((x * y) ** 2, 0), x >= 0)
    result = solve_semialgebraic(formula, (x, y), count=0)
    steps = result.diagnostics.get("planner_steps", ())
    ideal_steps = [step for step in steps if step["method"] == "equality_ideal_analysis"]
    assert ideal_steps
    assert ideal_steps[-1]["metadata"]["dimension"] == 1
    rur_steps = [step for step in steps if step["method"] == "rational_univariate"]
    assert rur_steps and rur_steps[-1]["accepted"] is False
