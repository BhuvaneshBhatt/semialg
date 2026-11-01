"""Presolve and lazy dispatch preserve polynomial syntax, roles and quantifiers."""

from unittest.mock import patch

import pytest
import sympy as sp

from semialg import quantifier_eliminate
from semialg.decision import equivalent
from semialg.errors import ExactComputationUnavailable
from semialg.formula import parse_formula
from semialg.partial.qe import lazy_find_inst_form, lazy_resolve_formula
from semialg.presolve import presolve_semialgebraic
from semialg.qe.complete import norm_internal_order, qe_by_complete_cad


@pytest.mark.parametrize("equality", [False, True])
@pytest.mark.parametrize("quantifier", ["exists", "forall"])
def test_nested_polynomial_disjunction(equality, quantifier):
    a, x, y = sp.symbols("a x y", real=True)
    outside = sp.Or(x < -1, x > 1)
    if equality:
        matrix = sp.And(a > 0, sp.Eq(x, y), outside.subs(x, y))
        prefix = ((quantifier, x), (quantifier, y))
    else:
        matrix = sp.And(a > 0, outside)
        prefix = ((quantifier, x),)
    result = quantifier_eliminate(matrix, prefix)
    assert equivalent(result, a > 0 if quantifier == "exists" else sp.false, (a,))
    presolved = presolve_semialgebraic(matrix, (a, x, y), eliminate=(x, y) if equality else ())
    assert not presolved.formula.has(sp.Abs)
    for atom in presolved.formula.atoms(sp.core.relational.Relational):
        assert sp.Poly(atom.lhs - atom.rhs, a, x, y).total_degree() <= 1


def test_presolved_quantifier_stays_bound():
    a, x, y = sp.symbols("a x y", real=True)
    matrix = sp.And(sp.Eq(x, y), sp.Eq(y * y, a))
    result = quantifier_eliminate(
        matrix, (("exists", x), ("exists", y)), variables=(a, x, y), return_result=True
    )
    assert equivalent(result.formula, a >= 0, (a,))
    assert result.free_variables == (a,)
    assert result.quantified_variables == (x, y)
    assert result.presolve_substitutions == ((x, y),)
    assert set(result.free_variables).isdisjoint(result.quantified_variables)


@pytest.mark.parametrize("kind,truth", [("forall", False), ("exists", True)])
def test_iterator_quantifiers(kind, truth):
    x = sp.Symbol("x", real=True)
    matrix = parse_formula(x > 0)
    prefix = ((kind, x),)
    resolved = lazy_resolve_formula((x,), iter(prefix), matrix)
    instance = lazy_find_inst_form((x,), matrix, quantifiers=iter(prefix))
    assert resolved.truth_value is truth
    assert instance.found is truth
    assert resolved.stats.quantifiers == prefix
    assert instance.stats.quantifiers == prefix
    if truth:
        assert instance.instance[x] > 0
    else:
        assert instance.instance is None
    assert norm_internal_order((x,), iter(prefix), matrix)[:3] == ((x,), (), (x,))


@pytest.mark.parametrize("strategy", ["auto", "planner", "partial", "lazy"])
def test_closed_qe_uses_partial_lifting(strategy):
    x, y = sp.symbols("x y", real=True)
    with patch(
        "semialg.qe.complete.decomp_collins_complete",
        side_effect=AssertionError("full CAD requested"),
    ):
        result = quantifier_eliminate(
            sp.Eq(y * y, x), (("forall", x), ("exists", y)), strategy=strategy, return_result=True
        )
    assert result.formula is sp.false
    assert result.method == "partial_cad_resolve"
    assert result.backend_result.backend == "partial-collins-lazy"
    assert result.backend_result.stats.stopped_early


def test_free_qe_keeps_complete_truth_coverage():
    a, x = sp.symbols("a x", real=True)
    result = quantifier_eliminate(
        sp.Eq(x * x, a), (("exists", x),), strategy="partial", return_result=True
    )
    assert equivalent(result.formula, a >= 0, (a,))
    assert result.free_variables == (a,)
    assert result.backend_result.truth_by_index == {(0,): False, (1,): True, (2,): True}


def test_lazy_failure_remains_unresolved():
    x, y = sp.symbols("x y", real=True)
    with patch(
        "semialg.partial.qe.lazy_resolve_formula",
        side_effect=ExactComputationUnavailable("forced lazy failure"),
    ):
        with pytest.raises(ExactComputationUnavailable, match="forced lazy failure"):
            quantifier_eliminate(sp.Eq(y * y, x), (("forall", x), ("exists", y)))


def test_complete_strategy_preserves_full_cad():
    x, y = sp.symbols("x y", real=True)
    result = quantifier_eliminate(
        sp.Eq(y * y, x), (("forall", x), ("exists", y)), strategy="collins", return_result=True
    )
    assert result.formula is sp.false
    assert result.method == "cad" and result.backend_result.cad.complete


@pytest.mark.slow
@pytest.mark.parametrize("first", ["exists", "forall"])
def test_four_level_prefix(first):
    from itertools import product

    variables = sp.symbols("w x y z", real=True)
    prefix = tuple(
        (first if i % 2 == 0 else "forall" if first == "exists" else "exists", variable)
        for i, variable in enumerate(variables)
    )
    predicate = sp.And(sp.Eq(variables[0], variables[1]), sp.Eq(variables[2], variables[3]))
    matrix = predicate
    for kind, variable in reversed(prefix):
        guard = sp.Eq(variable * (variable - 1), 0)
        matrix = sp.And(guard, matrix) if kind == "exists" else sp.Or(sp.Not(guard), matrix)
    outcomes = {
        values: bool(predicate.subs(dict(zip(variables, values, strict=True))))
        for values in product((0, 1), repeat=4)
    }
    for kind, _variable in reversed(prefix):
        keys = {key[:-1] for key in outcomes}
        outcomes = {
            key: (any if kind == "exists" else all)(outcomes[(*key, v)] for v in (0, 1))
            for key in keys
        }
    expected = outcomes[()]
    parsed = parse_formula(matrix)
    assert lazy_resolve_formula(variables, iter(prefix), parsed).truth_value is expected
    result = qe_by_complete_cad(variables, prefix, parsed, return_result=True)
    assert result.truth_value is expected
    assert bool(quantifier_eliminate(matrix, prefix)) is expected
