"""Common graph comparison preserves scopes, exact provenance and function obligations."""

from dataclasses import replace

import pytest
import sympy as sp

from semialg.function_graph import semialgebraic_function_graph
from semialg.function_verification import compare_real_graphs
from semialg.qe import quantifier_eliminate

x, y, a = sp.symbols("x y a", real=True)


def test_relation_and_function_obligations():
    graph = sp.Eq(y * y, x)
    result = compare_real_graphs(
        graph, graph, input=x, output=y, domain=x > 0, check_functions=True
    )
    assert result.equivalent and result.certified
    assert all(item.total and not item.single_valued for item in result.function_checks)
    assert result.function_verified is False


def test_auxiliary_scope():
    target = semialgebraic_function_graph(sp.sqrt(x), y)
    result = compare_real_graphs(
        target,
        sp.And(sp.Eq(y * y, x), y >= 0),
        input=x,
        output=y,
        domain=x >= 0,
        check_functions=True,
    )
    assert result.equivalent and result.certified and result.function_verified
    assert not set(target.auxiliary_variables) & result.target_formula.free_symbols


def test_projection_provenance():
    graph = quantifier_eliminate(sp.Eq(y, x), [], variables=(x, y), return_result=True)
    result = compare_real_graphs(graph, graph, input=x, output=y)
    assert result.eliminations[0] is graph and result.eliminations[1] is graph
    uncertain = replace(graph, certified=False)
    result = compare_real_graphs(uncertain, graph, input=x, output=y)
    assert result.equivalent and not result.certified


def test_counterexample():
    result = compare_real_graphs(
        sp.Eq(y, x), sp.Eq(y, x + 1), input=x, output=y, find_counterexample=True
    )
    assert not result.equivalent and result.certified
    point = result.comparison.counterexample
    assert point is not None
    assert sp.Xor(sp.Eq(y, x), sp.Eq(y, x + 1)).subs(point) is sp.true


def test_parameter_scope():
    result = compare_real_graphs(
        sp.Eq(a * y, x),
        sp.Eq(y, x / a),
        input=x,
        output=y,
        parameters=(a,),
        domain=sp.Ne(a, 0),
        check_functions=True,
    )
    assert result.equivalent and result.certified and result.function_verified
    projected = quantifier_eliminate(sp.Eq(y, a + x), [], variables=(a, x, y), return_result=True)
    with pytest.raises(ValueError, match="undeclared"):
        compare_real_graphs(projected, sp.Eq(y, x), input=x, output=y)
    with pytest.raises(ValueError, match="domain"):
        compare_real_graphs(sp.Eq(y, x), sp.Eq(y, x), input=x, output=y, domain=y > 0)
