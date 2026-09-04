import sympy as sp

from semialg import function_range, semialgebraic_minimize


def test_parametric_affine_range_eliminates_without_second_qe():
    x, a, t = sp.symbols("x a t", real=True)
    result = function_range(
        x,
        sp.And(x >= 0, x <= a),
        [x],
        parameters=[a],
        value_symbol=t,
        return_stratified=True,
        eliminate_quantifiers=True,
    )
    assert all(branch.value.quantifier_free for branch in result.branches)
    assert all(
        branch.value.method == "parametric_direct_range_relation" for branch in result.branches
    )
    positive = next(branch.value for branch in result.branches if branch.condition == (a > 0))
    assert (
        sp.simplify_logic(sp.Equivalent(positive.formula, sp.And(a > 0, t >= 0, t <= a))) is sp.true
    )


def test_parametric_affine_optimum_eliminates_without_second_qe():
    x, a = sp.symbols("x a", real=True)
    result = semialgebraic_minimize(
        x,
        sp.And(x >= 0, x <= a),
        [x],
        parameters=[a],
        return_stratified=True,
        eliminate_quantifiers=True,
        return_result=True,
    )
    assert all(branch.value.quantifier_free for branch in result.branches)
    assert all(
        branch.value.method == "parametric_direct_optimum_relation" for branch in result.branches
    )
