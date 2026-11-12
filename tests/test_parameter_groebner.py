"""Native generic bases must not remove exceptional real parameter fibers."""

import pytest
import sympy as sp
from algroots.native_f4 import compiled_available

from semialg import solvability_conditions
from semialg.algebraic.parameter_groebner import guarded_equality_preprocessing

x, y, t, u = sp.symbols("x y t u", real=True)


@pytest.mark.parametrize("backend", ["native"])
@pytest.mark.parametrize(
    "constraints,truth",
    [
        ([sp.Eq(t * x, 0), sp.Eq(x, 1)], [True, False, False]),
        ([sp.Eq(t * x - 1, 0)], [False, True, True]),
        ([sp.Eq(t * x, 0), sp.Gt(y, 0)], [True, True, True]),
    ],
)
def test_exceptional_solvability(backend, constraints, truth):
    if backend == "native" and not compiled_available():
        pytest.skip("native F4 is not installed")
    result = solvability_conditions(
        constraints, (x, y), (t,), groebner_backend=backend, return_result=True
    )
    assert result.diagnostics["exceptional_branch_preserved"]
    for value, expected in zip([0, 1, -1], truth, strict=True):
        assert bool(sp.simplify(result.formula.subs(t, value))) is expected


@pytest.mark.skipif(not compiled_available(), reason="native F4 extension unavailable")
def test_two_parameter_exception():
    result = solvability_conditions(
        [sp.Eq((t - u) * x, 0), sp.Eq(x, 1)], (x,), (t, u), groebner_backend="native"
    )
    assert result.subs({t: 1, u: 1}) == sp.true
    assert result.subs({t: 1, u: 0}) == sp.false


@pytest.mark.skipif(not compiled_available(), reason="native F4 extension unavailable")
def test_algebraic_parameter_exception():
    result = solvability_conditions(
        [sp.Eq((t - sp.sqrt(2)) * x, 0), sp.Eq(x, 1)], (x,), (t,), groebner_backend="native"
    )
    assert sp.simplify(result.subs(t, sp.sqrt(2))) == sp.true
    assert sp.simplify(result.subs(t, -sp.sqrt(2))) == sp.false


@pytest.mark.skipif(not compiled_available(), reason="native F4 extension unavailable")
def test_retained_inequalities():
    rewritten, diagnostics = guarded_equality_preprocessing(
        sp.And(sp.Eq(t * x, 0), x > 2), (x,), (t,), backend="native"
    )
    assert diagnostics["groebner_guard"] == t
    assert rewritten.subs({t: 0, x: 3}) == sp.true
    assert rewritten.subs({t: 0, x: 1}) == sp.false
    assert rewritten.subs({t: 1, x: 3}) == sp.false


@pytest.mark.parametrize("backend", ["sympy", "native"])
def test_visible_coefficient_pole(backend):
    if backend == "native" and not compiled_available():
        pytest.skip("native F4 is not installed")
    held = sp.Mul(x, sp.Mul(t, sp.Pow(t, -1, evaluate=False), evaluate=False), evaluate=False)
    result = solvability_conditions(
        sp.Eq(held, 0, evaluate=False), (x,), (t,), groebner_backend=backend
    )
    assert result.subs(t, 0) == sp.false
    assert result.subs(t, 1) == sp.true


@pytest.mark.parametrize("backend", ["sympy", "native"])
def test_rational_coordinate_domain(backend):
    if backend == "native" and not compiled_available():
        pytest.skip("native F4 is not installed")
    result = solvability_conditions(
        [sp.Eq(x / (x - 1), 0), sp.Eq(t * x, 0)], (x,), (t,), groebner_backend=backend
    )
    assert result == sp.true
    impossible = solvability_conditions(
        [sp.Eq(x / (x - 1), 0), sp.Eq(x, 1)], (x,), (t,), groebner_backend=backend
    )
    assert impossible == sp.false
