"""Empty witness searches are not certificates of infeasibility."""

from unittest.mock import patch

import pytest
import sympy as sp

from semialg import find_instance
from semialg.decomposition.components import component_instances
from semialg.errors import ExactComputationUnavailable
from semialg.instances.real_utils import FallbackInstanceResult


def test_inapplicable_ars():
    x = sp.Symbol("x", real=True)
    result = find_instance(x > 0, (x,), strategy="ars", return_result=True)
    assert result.status == "unknown"
    assert result.instances == ()
    assert bool((x > 0).subs(x, 1))


@pytest.mark.parametrize("strategy", ["heuristic", "random", "fallback"])
def test_exhausted_search(strategy):
    x = sp.Symbol("x", real=True)
    incomplete = FallbackInstanceResult((), "unknown", "forced_exhaustion", ())
    with patch("semialg.solve.find_instance.find_real_witnesses", return_value=incomplete):
        result = find_instance(x > 0, (x,), strategy=strategy, return_result=True)
    assert result.status == "unknown" and result.instances == ()


def test_cad_unknown_propagates():
    x = sp.Symbol("x", real=True)
    formula = sum(sp.sqrt(x + i) for i in range(4)) > 0
    assert formula.subs(x, 0) is sp.true
    components = component_instances(formula, (x,), return_result=True)
    assert components.status == "unknown" and components.instances == ()
    result = find_instance(formula, (x,), strategy="cad", return_result=True)
    assert result.status == "unknown" and result.instances == ()
    with pytest.raises(ExactComputationUnavailable):
        component_instances(formula, (x,))


def test_complex_tuple():
    x, y = sp.symbols("x y")
    condition = sp.And(sp.Eq(x, 1), sp.Eq(y, 2))
    result = find_instance(condition, (x, y), domain="complexes", return_result=True)
    assert result.status == "sat" and result.instances == ({x: 1, y: 2},)
    assert condition.subs(result.first()) is sp.true


def test_nonfinite_complex_set():
    z = sp.Symbol("z")
    result = find_instance(sp.Eq(sp.exp(z), 1), (z,), domain="complexes", return_result=True)
    assert result.status == "unknown" and result.instances == ()
    assert sp.Eq(sp.exp(z), 1).subs(z, 0) is sp.true


@pytest.mark.parametrize("domain", ["reals", "complexes", "integers", "booleans"])
def test_certified_empty(domain):
    x = sp.Symbol("x")
    result = find_instance(sp.false, (x,), domain=domain, return_result=True)
    assert result.status == "unsat" and result.instances == ()
