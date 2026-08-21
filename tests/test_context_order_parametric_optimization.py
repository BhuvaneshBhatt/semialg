from __future__ import annotations

import sympy as sp

import semialg
from semialg import (
    ExactComputationContext,
    ParameterStratifiedResult,
    ParametricFunctionRangeResult,
    ParametricOptimizationResult,
    computation_context,
    current_computation_context,
    function_range,
    semialgebraic_minimize,
)
from semialg.algebraic.roots import isolate_real_roots
from semialg.algebraic.signs import sign_at_sample
from semialg.optimization import _pruned_active_subsets
from semialg.planner.features import ProblemFeatures
from semialg.planner.heuristics import candidate_variable_orders, score_variable_order


def test_computation_context_reuses_algebraic_work_and_is_scoped() -> None:
    x = sp.Symbol("x", real=True)
    context = ExactComputationContext()
    assert current_computation_context() is None
    with computation_context(context):
        root = isolate_real_roots(x**2 - 2, x)[0]
        assert sign_at_sample(x + 3, (root,)) == 1
        assert sign_at_sample(x + 3, (root,)) == 1
        stats = context.stats()
        assert stats["algebraic.roots"]["size"] >= 1
        assert stats["algebraic.signs"]["hits"] >= 1
        assert current_computation_context() is context
    assert current_computation_context() is None


def test_nested_computation_context_reuses_outer_scope() -> None:
    outer = ExactComputationContext()
    inner = ExactComputationContext()
    with computation_context(outer):
        with computation_context(inner) as active:
            assert active is outer
            assert current_computation_context() is outer


def test_order_scoring_reports_arithmetic_complexity_and_pilot_lifting() -> None:
    x, y, z = sp.symbols("x y z", real=True)
    polys = (z**2 - x, y**2 + z - 1, x**2 - 2)
    features = ProblemFeatures(variables=(x, y, z), num_polynomials=len(polys))
    scores = candidate_variable_orders(features, polys, limit=6)
    measured = [score for score in scores if score.projection_poly_count is not None]
    assert measured
    assert all(score.estimated_alg_degree is not None for score in measured)
    assert all(score.coefficient_height_bits is not None for score in measured)
    piloted = [score for score in scores if score.pilot_lifting_roots is not None]
    assert piloted
    assert all(
        score.pilot_cell_count is not None and score.pilot_cell_count >= 1 for score in piloted
    )


def test_coefficient_height_contributes_to_order_score_metadata() -> None:
    x, y = sp.symbols("x y", real=True)
    small = score_variable_order(
        (x, y), (x + y + 1,), include_projection=True, include_lifting=True
    )
    large = score_variable_order(
        (x, y),
        (x + y + 2**80,),
        include_projection=True,
        include_lifting=True,
    )
    assert small.coefficient_height_bits is not None
    assert large.coefficient_height_bits is not None
    assert large.coefficient_height_bits > small.coefficient_height_bits


def test_active_set_semialgebraic_pruning_removes_infeasible_strict_boundary() -> None:
    x = sp.Symbol("x", real=True)
    subsets = _pruned_active_subsets((), (x,), (x,), x > 0)
    assert subsets == ((),)


def test_fully_eliminated_equality_reconstructs_optimizer_point() -> None:
    x = sp.Symbol("x", real=True)
    result = semialgebraic_minimize(x**2, [sp.Eq(x, 2)], [x])
    assert result.value == 4
    assert result.attained
    assert result.point == {x: 2}


def test_parametric_optimization_returns_first_class_exact_strata() -> None:
    x, a = sp.symbols("x a", real=True)
    result = semialgebraic_minimize(
        x**2,
        [sp.Eq(x, a)],
        [x],
        parameters=[a],
        return_stratified=True,
    )
    assert isinstance(result, ParameterStratifiedResult)
    assert result.parameters == (a,)
    branch = result.branches[0]
    assert isinstance(branch.value, ParametricOptimizationResult)
    assert branch.value.kind == "min"
    assert branch.value.parameters == (a,)
    assert branch.value.certified
    assert branch.value.quantifier_free is False
    assert [kind for kind, _ in branch.value.quantifiers] == ["forall", "forall", "exists"]
    assert a in branch.value.formula.free_symbols
    assert branch.value.value_symbol in branch.value.formula.free_symbols


def test_parametric_range_returns_first_class_exact_strata() -> None:
    x, a = sp.symbols("x a", real=True)
    result = function_range(
        x + a,
        [x >= 0, x <= 1],
        [x],
        parameters=[a],
        return_stratified=True,
    )
    assert isinstance(result, ParameterStratifiedResult)
    branch = result.branches[0]
    assert isinstance(branch.value, ParametricFunctionRangeResult)
    assert branch.value.parameters == (a,)
    assert branch.value.quantifiers[0][0] == "exists"
    quantified_x = branch.value.quantifiers[0][1]
    formula = branch.value.formula
    assert quantified_x in formula.free_symbols
    assert a in formula.free_symbols
    assert branch.value.value_symbol in formula.free_symbols


def test_parametric_string_parameter_preserves_existing_symbol_identity() -> None:
    x = sp.Symbol("x")
    a = sp.Symbol("a")
    result = function_range(
        x + a,
        [x >= 0, x <= 1],
        [x],
        parameters=["a"],
        return_stratified=True,
    )
    assert result.parameters == (a,)


def test_new_public_exports_resolve() -> None:
    for name in (
        "ExactComputationContext",
        "computation_context",
        "current_computation_context",
        "ParametricOptimizationResult",
        "ParametricFunctionRangeResult",
    ):
        assert getattr(semialg, name) is not None
