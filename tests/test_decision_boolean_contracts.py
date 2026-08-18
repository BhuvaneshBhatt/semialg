from __future__ import annotations

import sympy as sp

from semialg import equivalent, implies, is_satisfiable, is_tautology, sample_point, sign_vector
from semialg.decision import (
    EquivalenceResult,
    ImplicationResult,
    SatisfiabilityResult,
    TautologyResult,
)


def test_satisfiability_result_has_validated_witness():
    x = sp.symbols("x", real=True)
    result = is_satisfiable((x >= 0) & (x <= 1), [x], return_result=True)
    assert isinstance(result, SatisfiabilityResult)
    assert result.satisfiable
    assert result.witness is not None
    assert bool(((x >= 0) & (x <= 1)).subs(result.witness))
    assert is_satisfiable((x >= 0) & (x <= 1), [x]) is True


def test_tautology_result_counterexample():
    x = sp.symbols("x", real=True)
    result = is_tautology(x**2 > 0, [x], return_result=True)
    assert isinstance(result, TautologyResult)
    assert not result.tautology
    assert result.counterexample is not None
    assert result.counterexample[x] == 0


def test_implication_result_counterexample():
    x = sp.symbols("x", real=True)
    result = implies(x >= 0, x > 0, [x], return_result=True)
    assert isinstance(result, ImplicationResult)
    assert not result.valid
    assert result.counterexample is not None
    assert result.counterexample[x] == 0
    assert implies((x >= 0) & (x <= 1), x**2 <= 1, [x]) is True


def test_equivalence_result_failed_direction():
    x = sp.symbols("x", real=True)
    result = equivalent(x**2 < 1, (x >= -1) & (x <= 1), [x], return_result=True)
    assert isinstance(result, EquivalenceResult)
    assert not result.equivalent
    assert result.failed_direction in {"rhs_implies_lhs", "both"}
    assert result.counterexample is not None


def test_sample_point_is_validated_and_sign_vector_as_dict():
    x = sp.symbols("x", real=True)
    point = sample_point((x >= 0) & (x <= 1), [x])
    assert point is not None
    assert bool(((x >= 0) & (x <= 1)).subs(point))
    signs = sign_vector([x, x - 1, x**2 - 1], {x: 0}, as_dict=True)
    assert signs[x] == 0
    assert signs[x - 1] == -1
    assert signs[x**2 - 1] == -1
