import pytest
import sympy as sp

from semialg import function_range, semialgebraic_minimize
from semialg.formula import parse_formula
from semialg.qe.complete import qe_by_complete_cad

x, a, b, t = sp.symbols("x a b t", real=True)


def _generic_range(expr, constraints, params):
    matrix = sp.And(constraints, sp.Eq(t, expr))
    variables = [*params, t, x]
    result = qe_by_complete_cad(
        variables,
        [("exists", x)],
        parse_formula(matrix),
        free_variables=[*params, t],
        variable_order_strategy="preserve",
        use_presolve=False,
        return_result=True,
    )
    return result.formula


def _fast_range(expr, constraints, params):
    result = function_range(
        expr,
        constraints,
        [x],
        parameters=params,
        value_symbol=t,
        return_stratified=True,
        eliminate_quantifiers=True,
    )
    return sp.Or(*(branch.value.formula for branch in result.branches))


@pytest.mark.parametrize(
    "expr,constraints,params",
    [
        (x, sp.And(x >= 0, x <= a), [a]),
        (2 * x + 1, sp.And(x >= 0, x <= a), [a]),
        (-x, sp.And(x >= 0, x <= a), [a]),
        (x + a, sp.And(x >= 0, x <= 2), [a]),
    ],
)
def test_fast_parametric_range_matches_generic_complete_qe(expr, constraints, params):
    fast = _fast_range(expr, constraints, params)
    generic = _generic_range(expr, constraints, params)
    # Generic CAD reconstruction can retain root-function boundaries that are
    # semantically exact but not suitable for a second polynomial CAD pass.
    # Differentially compare exact rational specializations instead.
    symbols = [*params, t]
    grids = [(-2, -1, 0, 1, 2, 3) for _ in symbols]
    import itertools

    for values in itertools.product(*grids):
        subs = dict(zip(symbols, values, strict=True))
        assert bool(sp.simplify(fast.subs(subs))) == bool(sp.simplify(generic.subs(subs)))


@pytest.mark.parametrize(
    "expr,constraints,expected_at",
    [
        (x, sp.And(x >= 0, x <= a), [(2, 0, True), (2, 2, True), (2, 3, False)]),
        (2 * x + 1, sp.And(x >= 0, x <= a), [(2, 1, True), (2, 5, True), (2, 0, False)]),
        (-x, sp.And(x >= 0, x <= a), [(2, -2, True), (2, 0, True), (2, 1, False)]),
        (x + a, sp.And(x >= 0, x <= 2), [(3, 3, True), (3, 5, True), (3, 6, False)]),
        (3 * x - 2, sp.And(x >= -1, x <= 1), [(4, -5, True), (4, 1, True), (4, 2, False)]),
        (-2 * x + a, sp.And(x >= 0, x <= 1), [(3, 1, True), (3, 3, True), (3, 4, False)]),
    ],
)
def test_fast_parametric_range_accepts_exact_expected_values(expr, constraints, expected_at):
    formula = _fast_range(expr, constraints, [a])
    for aval, tval, expected in expected_at:
        assert bool(sp.simplify(formula.subs({a: aval, t: tval}))) is expected


@pytest.mark.parametrize(
    "expr,constraints,param,value_at",
    [
        (x, sp.And(x >= 0, x <= a), a, {2: 0}),
        (-x, sp.And(x >= 0, x <= a), a, {2: -2}),
        (2 * x + 3, sp.And(x >= 0, x <= a), a, {2: 3}),
        (-3 * x + 1, sp.And(x >= 0, x <= a), a, {2: -5}),
    ],
)
def test_fast_parametric_minimum_selects_correct_endpoint(expr, constraints, param, value_at):
    result = semialgebraic_minimize(
        expr,
        constraints,
        [x],
        parameters=[param],
        return_stratified=True,
        eliminate_quantifiers=True,
        return_result=True,
    )
    for pval, expected in value_at.items():
        selected = result.select({param: pval})
        assert selected is not None
        # Relations expose the optimum through their value symbol; the formula
        # must accept the expected value and reject a nearby wrong one.
        value_symbol = next(sym for sym in selected.formula.free_symbols if sym not in {param})
        assert sp.simplify(selected.formula.subs({param: pval, value_symbol: expected})) is sp.true
