import pickle

import sympy as sp

from semialg import is_satisfiable, solve_semialgebraic
from semialg.errors import UnsupportedFragmentError
from semialg.parameters import SolvabilityConditionsResult


def test_solution_result_formula_and_membership_contract():
    x = sp.Symbol("x", real=True)
    result = solve_semialgebraic(sp.And(x >= 0, x <= 2), [x])
    assert result.satisfiable is True
    formula = result.as_formula()
    assert is_satisfiable(formula, [x])
    if result.samples:
        assert result.contains(dict(result.samples[0]))


def test_parameter_result_dataclass_is_picklable_and_stable():
    a = sp.Symbol("a", real=True)
    result = SolvabilityConditionsResult(a >= 0, a >= 0, (), (a,))
    restored = pickle.loads(pickle.dumps(result))
    assert restored == result
    assert restored.parameters == (a,)


def test_public_unsupported_exception_is_distinct_from_bad_math_answer():
    exc = UnsupportedFragmentError("unsupported")
    assert isinstance(exc, ValueError)
    assert "unsupported" in str(exc)
