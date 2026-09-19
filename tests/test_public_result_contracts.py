import pickle

import sympy as sp

from semialg import is_satisfiable, solvability_conditions, solve_semialgebraic
from semialg.errors import UnsupportedFragmentError


def test_solution_result_formula_and_membership_contract():
    x = sp.Symbol("x", real=True)
    result = solve_semialgebraic(sp.And(x >= 0, x <= 2), [x])
    assert result.satisfiable is True
    formula = result.as_formula()
    assert is_satisfiable(formula, [x])
    if result.samples:
        assert result.contains(dict(result.samples[0]))


def test_parameter_result_from_public_operation_is_picklable_and_stable():
    x, a = sp.symbols("x a", real=True)
    result = solvability_conditions(sp.Eq(x, a), (x,), (a,), return_result=True)
    restored = pickle.loads(pickle.dumps(result))
    assert restored == result
    assert restored.parameters == (a,)
    assert restored.variables == (x,)
    assert restored.formula is sp.true


def test_public_unsupported_exception_is_distinct_from_bad_math_answer():
    exc = UnsupportedFragmentError("unsupported")
    assert isinstance(exc, ValueError)
    assert "unsupported" in str(exc)
