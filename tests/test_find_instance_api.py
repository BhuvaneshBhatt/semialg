import sympy as sp

from semialg import InstanceResult, find_instance, find_instance_text
from semialg.status import SolverStatus


def test_find_instance_01():
    x = sp.Symbol("x", real=True)
    result = find_instance(sp.Eq(x**2, 2), [x], count=2)
    assert isinstance(result, InstanceResult)
    assert result.status in {SolverStatus.SAT, SolverStatus.UNSAT}
    assert result.domain.value == "reals"


def test_find_instance_02():
    result = find_instance_text("x^2 = 2", variables=["x"], count=2)
    assert isinstance(result, InstanceResult)
    assert result.variables


def test_find_instance_03():
    x = sp.Symbol("x")
    result = find_instance(sp.Eq(x**2 + 1, 0), [x], domain="complexes", count=2)
    assert isinstance(result, InstanceResult)
    assert result.domain.value == "complexes"


def test_find_instance_04():
    x, y = sp.symbols("x y")
    result = find_instance(x | y, [x, y], domain="booleans", count=2)
    assert len(result.instances) == 2
