from __future__ import annotations

import pytest
import sympy as sp

from semialg import OptimizationResult, semialgebraic_maximize, semialgebraic_minimize


def test_exact_pipeline_interior_stationary_point() -> None:
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_minimize(x**2 + 2 * y**2, [x**2 + y**2 <= 4], [x, y])
    assert isinstance(result, OptimizationResult)
    assert result.value == 0
    assert result.point == {x: 0, y: 0}
    assert result.attained is True
    assert result.certified is True


def test_exact_pipeline_active_inequality_kkt() -> None:
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_minimize(x**2 + y**2, [x + y >= 1], [x, y])
    assert result.value == sp.Rational(1, 2)
    assert result.point == {x: sp.Rational(1, 2), y: sp.Rational(1, 2)}
    assert result.certified is True
    assert result.diagnostics["global_certificate"] == "complete_cad_no_better_point"


def test_exact_pipeline_multiple_active_boundaries_vertex() -> None:
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_minimize(x + 2 * y, [x >= 0, y >= 0, x <= 2, y <= 3], [x, y])
    assert result.value == 0
    assert result.point == {x: 0, y: 0}
    assert result.certified is True


def test_exact_pipeline_equality_constraint() -> None:
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_maximize(x, [sp.Eq(x**2 + y**2, 1)], [x, y])
    assert result.value == 1
    assert result.point == {x: 1, y: 0}
    assert result.certified is True


def test_exact_pipeline_singular_active_locus() -> None:
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_minimize(x, [sp.Eq(x**2 + y**2, 0)], [x, y])
    assert result.value == 0
    assert result.point == {x: 0, y: 0}
    assert result.certified is True


def test_exact_pipeline_open_boundary_reports_unattained_infimum() -> None:
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_minimize(x, [x**2 + y**2 < 1], [x, y])
    assert result.value == -1
    assert result.attained is False
    assert result.points == ()
    assert result.certified is True
    assert result.method.endswith("cad_range_certificate")


@pytest.mark.slow
def test_exact_pipeline_three_variables() -> None:
    x, y, z = sp.symbols("x y z", real=True)
    result = semialgebraic_maximize(x + y + z, [x**2 + y**2 + z**2 <= 1], [x, y, z])
    assert result.value == sp.sqrt(3)
    assert result.point == {x: sp.sqrt(3) / 3, y: sp.sqrt(3) / 3, z: sp.sqrt(3) / 3}
    assert result.certified is True
    assert result.method.endswith("cad_decision_certificate")


def test_exact_pipeline_disjunctive_domain() -> None:
    x = sp.symbols("x", real=True)
    domain = sp.Or(sp.And(x >= -3, x <= -2), sp.And(x >= 1, x <= 4))
    result = semialgebraic_minimize(x**2, domain, [x])
    assert result.value == 1
    assert result.point == {x: 1}
    assert result.attained is True
    assert result.certified is True


def test_exact_pipeline_algebraic_optimum_and_rur_candidates() -> None:
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_minimize(x, [sp.Eq(x**2, 2), sp.Eq(y, 0)], [x, y])
    assert result.value == -sp.sqrt(2)
    assert result.point == {x: -sp.sqrt(2), y: 0}
    assert result.certified is True


def test_exact_pipeline_return_value_api() -> None:
    x, y = sp.symbols("x y", real=True)
    value = semialgebraic_maximize(x * y, [x >= 0, y >= 0, x + y <= 1], [x, y], return_result=False)
    assert value == sp.Rational(1, 4)


def test_optimization_result_certified_field_is_public() -> None:
    x = sp.symbols("x", real=True)
    result = semialgebraic_minimize(x**2, [x >= 2], [x])
    assert result.certified is True
    assert isinstance(result.diagnostics, dict)


def test_string_variable_name_reuses_problem_symbol_without_assumption_duplicate() -> None:
    x = sp.Symbol("x")
    result = semialgebraic_minimize(x**2, [x >= 1], ["x"])
    assert result.variables == (x,)
    assert result.value == 1
    assert result.attained
    assert result.point is not None
    assert result.point[x] == 1


def test_disjunctive_optimization_requires_every_branch_to_be_certified(monkeypatch) -> None:
    import semialg.optimization as optimization_module

    x = sp.Symbol("x", real=True)
    calls = iter(
        (
            optimization_module.OptimizationResult(
                x, (x,), sp.Integer(0), ({x: 0},), True, "min", certified=True
            ),
            optimization_module.OptimizationResult(
                x, (x,), sp.Integer(2), ({x: 2},), True, "min", certified=False
            ),
        )
    )

    monkeypatch.setattr(
        optimization_module, "_optimize_conjunction", lambda *args, **kwargs: next(calls)
    )
    result = optimization_module.semialgebraic_minimize(
        x,
        sp.Or(sp.Eq(x, 0), sp.Eq(x, 2)),
        [x],
    )
    assert result.value == 0
    assert not result.certified
    assert result.diagnostics["all_branches_certified"] is False
