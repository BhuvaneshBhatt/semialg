from __future__ import annotations

import warnings

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

from semialg import BooleanSimplificationResult, equivalent, simplify_boole


def test_simplify_boole_canonicalizes_relational_atoms():
    x = sp.Symbol("x", real=True)

    assert simplify_boole(2 * x - 2 >= 0, [x]) == (x >= 1)
    assert simplify_boole(-x <= -1, [x]) == (x >= 1)


def test_simplify_boole_univariate_interval_normal_form():
    x = sp.Symbol("x", real=True)

    simplified = simplify_boole((x**2 <= 1) & (x >= 0), [x])
    assert equivalent(simplified, (x >= 0) & (x <= 1), [x])
    assert simplified == ((x >= 0) & (x <= 1))


def test_simplify_boole_univariate_interval_tautology():
    x = sp.Symbol("x", real=True)

    assert simplify_boole((x**2 > 1) | ((x >= -1) & (x <= 1)), [x]) == sp.true


def test_simplify_boole_absorbs_semantically_redundant_disjunct():
    x = sp.Symbol("x", real=True)

    assert simplify_boole((x >= 0) | (x > 1), [x]) == (x >= 0)


def test_simplify_boole_return_result_contract():
    x = sp.Symbol("x", real=True)

    result = simplify_boole((x >= 0) & (x > 1), [x], return_result=True)
    assert isinstance(result, BooleanSimplificationResult)
    assert result.formula == (x > 1)
    assert result.variables == (x,)
    assert result.original == ((x >= 0) & (x > 1))


def test_simplify_boole_does_not_call_scalar_simplify_on_boolean_formula():
    x = sp.Symbol("x", real=True)
    formula = (x >= 0) & (x <= 1)

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        simplify_boole(formula, [x])

    assert not [w for w in captured if issubclass(w.category, SymPyDeprecationWarning)]
