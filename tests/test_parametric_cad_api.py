from __future__ import annotations

import pytest
import sympy as sp

from semialg.decomposition import ParametricCADFunction, parametric_cad, parametric_cad_text


def test_parametric_cad_01():
    a, x = sp.symbols("a x", real=True)
    result = parametric_cad(sp.Eq(a * x - 1, 0), variables=[x], parameters=[a])
    assert result.status == "complete"
    assert bool(result.generic_formula.subs({a: 2}))
    assert not bool(result.generic_formula.subs({a: 0}))
    assert bool(result.exceptional_formula.subs({a: 0}))
    assert result.generic_cases
    assert result.all_generic_cases
    assert result.exceptional_cases
    assert all(case.param_condition is not sp.false for case in result.cases)


@pytest.mark.slow
def test_parametric_cad_02():
    cases = parametric_cad_text(
        "a*x - 1 == 0", variables=["x"], parameters=["a"], output="cases", return_result=False
    )
    assert cases
    assert any(case.exceptional for case in cases)


@pytest.mark.slow
def test_parametric_cad_03():
    a, x = sp.symbols("a x", real=True)
    fn = parametric_cad(sp.Eq(a * x - 1, 0), variables=[x], parameters=[a], output="function")
    assert isinstance(fn, ParametricCADFunction)
    assert bool(fn({a: 2}).subs({x: sp.Rational(1, 2)}))
    assert fn.exceptional({a: 0})


def test_parametric_cad_04():
    assert callable(parametric_cad)


def test_parametric_cad_reuses_cases():
    a, x = sp.symbols("a x", real=True)
    fn = parametric_cad(sp.Eq(x**2, a), variables=[x], parameters=[a], output="function")
    fiber = fn.specialize({a: sp.Integer(4)})
    assert bool(fiber.subs({x: 2}))
    assert not bool(fiber.subs({x: 3}))
    normal = fn.normal_form()
    assert bool(normal.subs({a: 4, x: -2}))
    assert not bool(normal.subs({a: -1, x: 0}))
