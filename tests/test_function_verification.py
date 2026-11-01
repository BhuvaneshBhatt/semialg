"""Behavioral contracts for exact real graph verification."""

import pytest
import sympy as sp

from semialg import equivalent, semialgebraic_function_graph
from semialg.equivalence import sym_diff_empty
from semialg.function_verification import verify_real_composition, verify_real_function

x, y, z, a, b, c, t = sp.symbols("x y z a b c t", real=True)


def test_fourth_root_graph_composition_is_exact():
    target = sp.And(sp.Eq(x**4, z), x >= 0)
    outer = sp.And(sp.Eq(x * x, y), x >= 0)
    inner = sp.And(sp.Eq(y * y, z), y >= 0)
    result = verify_real_composition(
        target, outer, inner, input=z, output=x, intermediate=y, domain=z >= 0
    )
    assert result.equivalent and result.certified
    assert all(r.certified for r in result.eliminations)
    assert equivalent(result.composition_formula, target, (x, z))


def test_missing_sign_selector_changes_the_real_branch():
    target = sp.And(sp.Eq(x**4, z), x >= 0)
    result = verify_real_composition(
        target, sp.Eq(x * x, y), sp.And(sp.Eq(y * y, z), y >= 0), input=z, output=x, intermediate=y
    )
    assert not result.equivalent
    assert result.comparison.counterexample is not None


def test_graph_objects_preserve_radicals_and_bind_auxiliaries():
    target = semialgebraic_function_graph(z ** sp.Rational(1, 4), x)
    outer = semialgebraic_function_graph(sp.sqrt(y), x)
    inner = semialgebraic_function_graph(sp.sqrt(z), y)
    result = verify_real_composition(target, outer, inner, input=z, output=x, intermediate=y)
    assert result.equivalent
    assert (
        not (
            set(target.auxiliary_variables)
            | set(outer.auxiliary_variables)
            | set(inner.auxiliary_variables)
        )
        & result.composition_formula.free_symbols
    )


def test_poles_are_preserved_when_comparing_rational_graphs():
    target = semialgebraic_function_graph(z, x)
    outer = semialgebraic_function_graph(1 / y, x)
    inner = semialgebraic_function_graph(1 / z, y)
    global_result = verify_real_composition(target, outer, inner, input=z, output=x, intermediate=y)
    assert not global_result.equivalent
    restricted = verify_real_composition(
        target, outer, inner, input=z, output=x, intermediate=y, domain=sp.Ne(z, 0)
    )
    assert restricted.equivalent


def test_parameter_domain_and_independent_auxiliary_scopes():
    target = semialgebraic_function_graph(z + a, x)
    outer = semialgebraic_function_graph(y + a, x)
    inner = semialgebraic_function_graph(z, y)
    result = verify_real_composition(
        target, outer, inner, input=z, output=x, intermediate=y, parameters=(a,), domain=a >= 0
    )
    assert result.equivalent


def test_graph_contract_rejects_undeclared_variables_and_capture():
    with pytest.raises(ValueError, match="undeclared"):
        verify_real_composition(
            sp.Eq(x, z), sp.Eq(x, y + z), sp.Eq(y, z), input=z, output=x, intermediate=y
        )
    with pytest.raises(ValueError, match="domain"):
        verify_real_composition(
            sp.Eq(x, z), sp.Eq(x, y), sp.Eq(y, z), input=z, output=x, intermediate=y, domain=x > 0
        )
    with pytest.raises(ValueError, match="target"):
        verify_real_composition(
            semialgebraic_function_graph(z, t),
            sp.Eq(x, y),
            sp.Eq(y, z),
            input=z,
            output=x,
            intermediate=y,
        )


def test_historical_equivalence_helper_cannot_prove_equivalence_by_sampling():
    result = sym_diff_empty(x < 4, sp.true, (x,))
    assert not result.equivalent and result.proof_attempted
    equal = sym_diff_empty(x * x <= 1, sp.And(x >= -1, x <= 1), (x,))
    assert equal.equivalent and equal.proof_attempted
    assert equal.method != "grid"


def test_real_graph_totality_uniqueness_and_domain_restriction():
    sqrt_graph = sp.And(sp.Eq(y * y, x), y >= 0)
    principal = verify_real_function(sqrt_graph, input=x, output=y, domain=x >= 0)
    assert principal.total and principal.single_valued and principal.certified
    partial = verify_real_function(sqrt_graph, input=x, output=y)
    assert not partial.total and partial.single_valued and partial.certified
    multivalued = verify_real_function(sp.Eq(y * y, x), input=x, output=y, domain=x > 0)
    assert multivalued.total and not multivalued.single_valued
    empty = verify_real_function(sp.false, input=x, output=y)
    assert not empty.total and empty.single_valued
    vacuous = verify_real_function(sp.false, input=x, output=y, domain=sp.false)
    assert vacuous.total and vacuous.single_valued


def test_real_rational_pole_and_parameter_totality():
    graph = sp.Eq(x * y, 1)
    unrestricted = verify_real_function(graph, input=x, output=y)
    assert not unrestricted.total and unrestricted.single_valued
    restricted = verify_real_function(graph, input=x, output=y, domain=sp.Ne(x, 0))
    assert restricted.total and restricted.single_valued
    parametric = verify_real_function(
        sp.Eq(a * y, x), input=x, output=y, parameters=(a,), domain=sp.Ne(a, 0)
    )
    assert parametric.total and parametric.single_valued


def test_graph_contract_validation():
    with pytest.raises(ValueError, match="domain"):
        verify_real_function(sp.Eq(y, x), input=x, output=y, domain=y > 0)
    with pytest.raises(ValueError, match="undeclared"):
        verify_real_function(sp.Eq(y, x + a), input=x, output=y)


def test_graph_object_auxiliaries_are_projected_before_function_checks():
    from semialg.function_graph import semialgebraic_function_graph

    graph = semialgebraic_function_graph(sp.sqrt(x), y)
    result = verify_real_function(graph, input=x, output=y, domain=x >= 0)
    assert result.total and result.single_valued and result.certified
    assert not (set(graph.auxiliary_variables) & result.formula.free_symbols)
