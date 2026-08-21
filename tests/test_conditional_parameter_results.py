from __future__ import annotations

import pytest
import sympy as sp

from semialg import (
    ConditionalBranch,
    ParameterStratifiedResult,
    conditional_result,
    parameterized_cylindrical_decomposition,
    root_count_conditions,
    solvability_conditions,
    verify_parameter_stratification,
)


def test_conditional_branch_applies_and_specializes() -> None:
    a = sp.symbols("a", real=True)
    branch = ConditionalBranch(a >= 0, a**2)
    assert branch.applies({a: 2}) is True
    assert branch.applies({a: -1}) is False
    specialized = branch.specialize({a: 3})
    assert specialized.condition is sp.true or specialized.condition == sp.true
    assert specialized.value == 9


def test_conditional_result_selects_disjoint_branch() -> None:
    a = sp.symbols("a", real=True)
    result = conditional_result([a], [(a < 0, -a), (a >= 0, a)])
    assert isinstance(result, ParameterStratifiedResult)
    assert result.select({a: -3}) == 3
    assert result.select({a: 4}) == 4


def test_conditional_result_requires_sufficient_assignments() -> None:
    a, b = sp.symbols("a b", real=True)
    result = conditional_result([a, b], [(a < b, 1), (a >= b, 2)])
    with pytest.raises(ValueError, match="insufficient"):
        result.select({a: 0})


def test_normalize_merges_equal_values() -> None:
    a = sp.symbols("a", real=True)
    result = conditional_result([a], [(a < 0, 7), (sp.Eq(a, 0), 7), (a > 0, 9)])
    assert result.stratum_count == 2
    assert sp.simplify_logic(result.condition_for_value(7) ^ (a <= 0)) is sp.false


def test_as_piecewise_and_partial_specialization() -> None:
    a, b = sp.symbols("a b", real=True)
    result = conditional_result([a, b], [(a < 0, b), (a >= 0, b + 1)])
    piecewise = result.as_piecewise()
    assert sp.simplify(piecewise.subs({a: -1, b: 3})) == 3
    partial = result.specialize({b: 5})
    assert partial.parameters == (a,)
    assert partial.select({a: 1}) == 6


def test_verify_parameter_stratification_accepts_partition() -> None:
    a = sp.symbols("a", real=True)
    result = conditional_result([a], [(a < 0, -1), (sp.Eq(a, 0), 0), (a > 0, 1)])
    cert = verify_parameter_stratification(result)
    assert cert.verify()
    assert cert.pairwise_disjoint
    assert cert.coverage_verified


def test_verify_parameter_stratification_detects_overlap() -> None:
    a = sp.symbols("a", real=True)
    result = conditional_result([a], [(a <= 1, 0), (a >= 0, 1)], normalize=False)
    cert = verify_parameter_stratification(result)
    assert not cert.verify()
    assert not cert.pairwise_disjoint
    assert cert.overlap_conditions


def test_solvability_conditions_first_class_strata() -> None:
    x, a = sp.symbols("x a", real=True)
    result = solvability_conditions(sp.Eq(x**2, a), [x], [a], return_stratified=True)
    assert isinstance(result, ParameterStratifiedResult)
    assert result.select({a: 2}) is True
    assert result.select({a: -2}) is False
    assert verify_parameter_stratification(result).verify()


def test_root_count_conditions_first_class_strata() -> None:
    x, a = sp.symbols("x a", real=True)
    result = root_count_conditions(x**2 + a, x, [a], return_stratified=True)
    assert isinstance(result, ParameterStratifiedResult)
    assert result.select({a: -1}) == 2
    assert result.select({a: 0}) == 1
    assert result.select({a: 1}) == 0
    assert verify_parameter_stratification(result).verify()


def test_result_objects_convert_to_stratified_results() -> None:
    x, a = sp.symbols("x a", real=True)
    solvability = solvability_conditions(sp.Eq(x**2, a), [x], [a], return_result=True)
    root_counts = root_count_conditions(x**2 + a, x, [a], return_result=True)
    assert solvability.as_stratified_result().select({a: 4}) is True
    assert root_counts.as_stratified_result().select({a: 4}) == 0


def test_parameterized_decomposition_exposes_guarded_strata_not_sample_fibers() -> None:
    x, a = sp.symbols("x a", real=True)
    decomposition = parameterized_cylindrical_decomposition(sp.Eq(x**2, a), [x], [a])
    result = decomposition.as_stratified_result()
    assert isinstance(result, ParameterStratifiedResult)
    selected = result.select({a: 1})
    assert selected.condition.subs(a, 1) is sp.true or selected.condition.subs(a, 1) == sp.true
    # The value is the ParameterStratum itself; representative ``solution`` is
    # explicitly not promoted to a symbolic fiber valid throughout the cell.
    assert hasattr(selected, "sample")


def test_public_conditional_exports_resolve() -> None:
    import semialg

    for name in (
        "ConditionalBranch",
        "ParameterStratifiedResult",
        "ParameterStratificationCertificate",
        "conditional_result",
        "verify_parameter_stratification",
    ):
        assert getattr(semialg, name) is not None


def test_string_assignments_resolve_to_original_parameter_symbol_assumptions() -> None:
    a = sp.Symbol("a")
    result = conditional_result([a], [(a < 0, -a), (a >= 0, a)])
    assert result.select({"a": -3}) == 3
    partial = result.specialize({"a": 4})
    assert partial.parameters == ()
    assert partial.select({}) == 4


def test_string_parameter_names_bind_to_symbols_already_used_by_branches() -> None:
    a = sp.Symbol("a")
    result = conditional_result(["a"], [(a < 0, -1), (a >= 0, 1)])
    assert result.parameters == (a,)
    assert result.select({"a": 2}) == 1
