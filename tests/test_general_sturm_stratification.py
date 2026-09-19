import sympy as sp
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg.parameters import root_count_conditions


def _distinct_real_root_count(expr: sp.Expr, variable: sp.Symbol) -> int:
    return len(tuple(dict.fromkeys(sp.real_roots(sp.Poly(expr, variable)))))


def test_degree_five_family_gets_complete_sturm_stratification() -> None:
    x, a = sp.symbols("x a", real=True)
    result = root_count_conditions(x**5 + a * x + 1, x, [a], return_result=True)

    assert result.method == "subresultant_sturm_parameter_stratification"
    assert result.classification.diagnostics["complete"] is True
    assert sp.Integer(-1) not in result.conditions_by_count
    assert set(result.conditions_by_count) == {sp.Integer(1), sp.Integer(2), sp.Integer(3)}


def test_degree_six_family_handles_repeated_root_stratum() -> None:
    x, a = sp.symbols("x a", real=True)
    family = (x**2 - a) ** 2 * (x**2 + 1)
    result = root_count_conditions(family, x, [a], return_result=True)

    assert result.classification.diagnostics["complete"] is True
    assert sp.Integer(-1) not in result.conditions_by_count
    assert result.condition_for_count(2).subs(a, 1) is sp.true
    assert result.condition_for_count(1).subs(a, 0) is sp.true
    assert result.condition_for_count(0).subs(a, -1) is sp.true


@given(st.integers(min_value=-8, max_value=8))
@settings(max_examples=17, deadline=None)
def test_degree_five_strata_agree_with_exact_specialization(a_value: int) -> None:
    x, a = sp.symbols("x a", real=True)
    family = x**5 + a * x + 1
    result = root_count_conditions(family, x, [a], return_result=True)
    expected = _distinct_real_root_count(family.subs(a, a_value), x)

    selected = [
        int(count)
        for count, condition in result.conditions_by_count.items()
        if sp.simplify(condition.subs(a, a_value)) is sp.true
    ]
    assert selected == [expected]
