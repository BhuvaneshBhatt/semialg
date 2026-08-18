from __future__ import annotations

import pytest
import sympy as sp

from semialg import GenericCADFunction, generic_cad, generic_cad_text


def test_generic_cad_01():
    a, x = sp.symbols("a x", real=True)
    result = generic_cad(sp.Eq(a * x - 1, 0), variables=[x], parameters=[a])
    assert result.status == "complete"
    assert bool(result.generic_formula.subs({a: 2}))
    assert not bool(result.generic_formula.subs({a: 0}))
    assert bool(result.exceptional_formula.subs({a: 0}))
    assert result.generic_cases
    assert result.all_generic_cases
    assert result.exceptional_cases
    assert all(case.param_condition is not sp.false for case in result.cases)


@pytest.mark.slow
def test_generic_cad_02():
    cases = generic_cad_text(
        "a*x - 1 == 0", variables=["x"], parameters=["a"], output="cases", return_result=False
    )
    assert cases
    assert any(case.exceptional for case in cases)


@pytest.mark.slow
def test_generic_cad_03():
    a, x = sp.symbols("a x", real=True)
    fn = generic_cad(
        sp.Eq(a * x - 1, 0), variables=[x], parameters=[a], output="function", return_result=False
    )
    assert isinstance(fn, GenericCADFunction)
    assert bool(fn({a: 2}).subs({x: sp.Rational(1, 2)}))
    assert fn.exceptional({a: 0})


def test_generic_cad_04():
    import semialg

    assert semialg.generic_cad is generic_cad
