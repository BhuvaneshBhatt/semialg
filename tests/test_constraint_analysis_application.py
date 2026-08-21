import pytest
import sympy as sp

from semialg.applications import analyze_constraint_redundancy, diagnose_feasible_set


def test_constraint_redundancy_finds_implied_constraint():
    x = sp.Symbol("x", real=True)
    result = analyze_constraint_redundancy([x >= 2, x >= 1, x <= 5], [x])
    assert result.certified
    assert result.feasible
    assert result.redundant_indices == (1,)
    assert result.essential_indices == (0, 2)
    assert result.witnesses[0] is not None


def test_constraint_redundancy_handles_no_redundancies():
    x = sp.Symbol("x", real=True)
    result = analyze_constraint_redundancy([x >= 0, x <= 1], [x])
    assert result.redundant_indices == ()
    assert result.essential_indices == (0, 1)


def test_feasible_set_diagnostics_returns_witness_and_redundancy():
    x = sp.Symbol("x", real=True)
    result = diagnose_feasible_set([x >= 2, x >= 1, x <= 5], [x])
    assert result.feasible
    assert result.witness is not None
    assert result.redundant_indices == (1,)
    assert result.conflict_indices == ()


def test_feasible_set_diagnostics_returns_irreducible_conflict():
    x = sp.Symbol("x", real=True)
    constraints = [x > 3, x < 2, x**2 <= 100]
    result = diagnose_feasible_set(constraints, [x])
    assert not result.feasible
    assert set(result.conflict_indices) == {0, 1}
    for index in result.conflict_indices:
        reduced = [
            constraint
            for pos, constraint in enumerate(constraints)
            if pos in result.conflict_indices and pos != index
        ]
        assert diagnose_feasible_set(reduced, [x], find_conflict=False).feasible


def test_constraint_analysis_rejects_empty_constraint_list():
    x = sp.Symbol("x", real=True)
    with pytest.raises(ValueError, match="at least one constraint"):
        analyze_constraint_redundancy([], [x])
    with pytest.raises(ValueError, match="at least one constraint"):
        diagnose_feasible_set([], [x])
