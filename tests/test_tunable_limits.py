from __future__ import annotations

import inspect

import pytest
import sympy as sp

from semialg import cad, sample_points, semialgebraic_minimize
from semialg.algebraic import configure_algebraic_cache_limits
from semialg.algebraic.cache import CACHE as ALGEBRAIC_CACHE
from semialg.cad_algorithms import configure_cad_cache_limits
from semialg.cad_algorithms.performance_cache import (
    COMPLETE_CADS,
    PROJECTION_STEPS,
    PROJECTION_TOWERS,
    SQUAREFREE_BASES,
)
from semialg.decomposition import CellSet, cad_text
from semialg.errors import ResourceLimitError
from semialg.planner.heuristics import candidate_variable_orders
from semialg.simplify import simplify_semialgebraic_formula


def test_cad_preprocess_aux_limit_never_turns_unknown_into_true():
    x = sp.Symbol("x", real=True)
    formula = sp.Abs(x) <= 1

    result = cad(formula, [x], max_preprocess_aux_vars=0, return_result=True)

    assert result.status == "unknown"
    assert result.formula == formula
    assert result.diagnostics["preprocess_elimination_limited"] is True
    assert result.diagnostics["max_preprocess_aux_vars"] == 0
    assert result.diagnostics["qe_formula"] is None
    with pytest.raises(ResourceLimitError):
        result.as_function()


def test_cad_preprocess_aux_limit_raises_when_structured_status_is_unavailable():
    x = sp.Symbol("x", real=True)
    formula = sp.Abs(x) <= 1

    with pytest.raises(ResourceLimitError):
        cad(formula, [x], max_preprocess_aux_vars=0, return_result=False)
    with pytest.raises(ResourceLimitError):
        cad(formula, [x], max_preprocess_aux_vars=0, strict=True)


def test_sampling_default_radius_and_random_denominator_are_user_tunable():
    x = sp.Symbol("x", real=True)

    assert (
        sample_points(
            sp.Eq(x, 2), [x], strategy="grid", grid_resolution=5, default_sampling_radius=1
        )
        == ()
    )
    assert sample_points(
        sp.Eq(x, 2), [x], strategy="grid", grid_resolution=5, default_sampling_radius=2
    ) == ({x: sp.Integer(2)},)

    point = sample_points(
        sp.true,
        [x],
        strategy="random",
        count=1,
        seed=7,
        random_attempts=1,
        max_random_denominator=7,
    )[0]
    assert sp.denom(point[x]) <= 7


def test_optimization_boolean_branch_limit_is_user_tunable():
    x = sp.Symbol("x", real=True)
    condition = sp.Or(x <= -1, x >= 1)

    with pytest.raises(NotImplementedError, match="max_boolean_branches=1"):
        semialgebraic_minimize(x**2, condition, [x], max_boolean_branches=1, return_result=True)


def test_simplification_and_planner_resource_guards_are_exposed():
    simplify_params = inspect.signature(simplify_semialgebraic_formula).parameters
    assert {
        "max_dnf_branches",
        "max_implication_atoms",
        "max_implication_vars",
        "max_implication_degree",
    } <= set(simplify_params)

    planner_params = inspect.signature(candidate_variable_orders).parameters
    assert {
        "limit",
        "exhaustive_var_limit",
        "projection_var_limit",
        "projection_shortlist",
    } <= set(planner_params)


def test_process_local_algebraic_and_cad_cache_capacities_are_configurable():
    old_alg = configure_algebraic_cache_limits()
    old_cad = configure_cad_cache_limits()
    try:
        new_alg = configure_algebraic_cache_limits(
            roots=17, signs=19, comparisons=23, specializations=29, rur=31
        )
        assert (new_alg.roots, new_alg.signs, new_alg.comparisons) == (17, 19, 23)
        assert ALGEBRAIC_CACHE.roots.maxsize == 17
        assert ALGEBRAIC_CACHE.specializations.maxsize == 29
        assert ALGEBRAIC_CACHE.rur.maxsize == 31

        new_cad = configure_cad_cache_limits(
            projection_towers=11,
            squarefree_bases=13,
            projection_steps=17,
            complete_cads=19,
        )
        assert new_cad.complete_cads == 19
        assert PROJECTION_TOWERS.maxsize == 11
        assert SQUAREFREE_BASES.maxsize == 13
        assert PROJECTION_STEPS.maxsize == 17
        assert COMPLETE_CADS.maxsize == 19
    finally:
        configure_algebraic_cache_limits(
            roots=old_alg.roots,
            signs=old_alg.signs,
            comparisons=old_alg.comparisons,
            specializations=old_alg.specializations,
            rur=old_alg.rur,
        )
        configure_cad_cache_limits(
            projection_towers=old_cad.projection_towers,
            squarefree_bases=old_cad.squarefree_bases,
            projection_steps=old_cad.projection_steps,
            complete_cads=old_cad.complete_cads,
        )


def test_cad_default_returns_requested_direct_representation():
    x = sp.Symbol("x", real=True)

    formula = cad(x**2 <= 1, [x])
    cells = cad(x**2 <= 1, [x], output="cells")

    assert isinstance(formula, sp.Basic)
    assert not hasattr(formula, "status")
    assert isinstance(cells, CellSet)

    text_formula = cad_text("x^2 <= 1", variables=["x"])
    text_cells = cad_text("x^2 <= 1", variables=["x"], output="cells")
    assert isinstance(text_formula, sp.Basic)
    assert isinstance(text_cells, CellSet)


def test_cad_unsupported_domain_cannot_be_mistaken_for_false_direct_output():
    x = sp.Symbol("x", real=True)

    with pytest.raises(NotImplementedError, match="only the real domain"):
        cad(x >= 0, [x], domain="integers")

    result = cad(x >= 0, [x], domain="integers", return_result=True)
    assert result.status == "unknown"
    assert result.formula == (x >= 0)
    assert "only the real domain" in result.diagnostics["reason"]


def test_cad_rejects_invalid_control_options_before_computation():
    x = sp.Symbol("x", real=True)

    with pytest.raises(ValueError, match="unsupported CAD output"):
        cad(x >= 0, [x], output="bogus")
    with pytest.raises(ValueError, match="unsupported CAD operation"):
        cad(x >= 0, [x], operation="bogus")
    with pytest.raises(ValueError, match="formula_form"):
        cad(x >= 0, [x], formula_form="bogus")
    with pytest.raises(ValueError, match="max_formula_terms"):
        cad(x >= 0, [x], max_formula_terms=0)
    with pytest.raises(ValueError, match="max_cells"):
        cad(x >= 0, [x], max_cells=0)
    with pytest.raises(ValueError, match="timeout"):
        cad(x >= 0, [x], timeout=0)
