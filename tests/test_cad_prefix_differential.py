"""Finite-domain oracles for short quantified prefixes and backend restarts."""

from unittest.mock import patch

import pytest
import sympy as sp

from semialg.errors import ExactComputationUnavailable, ProjectionNullification
from semialg.formula import parse_formula
from semialg.partial.qe import lazy_resolve_formula
from semialg.qe.complete import qe_by_complete_cad


@pytest.mark.parametrize("length", [1, 2, 3])
@pytest.mark.parametrize("first", ["exists", "forall"])
def test_alternating_prefix(length, first):
    variables = sp.symbols("x:z", real=True)[:length]
    quantifiers = tuple(
        ((first if i % 2 == 0 else "forall" if first == "exists" else "exists"), v)
        for i, v in enumerate(variables)
    )
    polynomial = sum((-1) ** i * v**2 for i, v in enumerate(variables))
    predicate = sp.Eq(polynomial, 0)
    matrix = predicate
    for quantifier, variable in reversed(quantifiers):
        guard = sp.Eq(variable * (variable - 1), 0)
        matrix = sp.And(guard, matrix) if quantifier == "exists" else sp.Or(sp.Not(guard), matrix)

    def finite_truth(position, assignment):
        if position == length:
            return bool(predicate.subs(assignment))
        quantifier, variable = quantifiers[position]
        values = [finite_truth(position + 1, {**assignment, variable: value}) for value in (0, 1)]
        return any(values) if quantifier == "exists" else all(values)

    expected = finite_truth(0, {})
    parsed = parse_formula(matrix)
    lazy = lazy_resolve_formula(variables, quantifiers, parsed)
    complete = qe_by_complete_cad(variables, quantifiers, parsed, return_result=True)
    assert lazy.truth_value is expected
    assert complete.truth_value is expected


@pytest.mark.parametrize("first", ["exists", "forall"])
def test_restart_retains_quantifiers(first):
    from semialg.partial.qe import _stack_for_parent

    x, y = sp.symbols("x y", real=True)
    quantifiers = ((first, x), ("forall" if first == "exists" else "exists", y))
    parsed = parse_formula(sp.Eq(y**2, x))
    calls = []

    def first_stack_fails(state, parent, level, **kwargs):
        calls.append(state.tower.metadata["projection"])
        if len(calls) == 1:
            raise ProjectionNullification("forced curtain")
        return _stack_for_parent(state, parent, level, **kwargs)

    with patch("semialg.partial.qe._stack_for_parent", side_effect=first_stack_fails):
        result = lazy_resolve_formula((x, y), quantifiers, parsed)
    # Neither one y works for all y, nor can negative x be a real square.
    assert result.truth_value is False
    assert result.backend == "partial-lazard-lazy"
    assert calls[0] != "lazard-complete"
    assert set(calls[1:]) == {"lazard-complete"}
    assert result.stats.quantifiers == quantifiers


def test_restart_failure_is_not_a_truth_value():
    x, y = sp.symbols("x y", real=True)
    with (
        patch(
            "semialg.partial.qe._stack_for_parent",
            side_effect=ProjectionNullification("forced curtain"),
        ),
        patch(
            "semialg.cad_algorithms.projection.lazard_complete.build_lazard_complete_projection",
            side_effect=ExactComputationUnavailable("forced projection failure"),
        ),
    ):
        with pytest.raises(ExactComputationUnavailable, match="forced projection failure"):
            lazy_resolve_formula(
                (x, y), (("forall", x), ("exists", y)), parse_formula(sp.Eq(y**2, x))
            )
