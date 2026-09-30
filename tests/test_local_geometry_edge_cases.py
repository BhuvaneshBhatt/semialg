import pytest
import sympy as sp

from semialg import (
    find_instance,
    local_components,
    local_germ,
    local_preimage,
    local_sign_strata,
    parameter_strata,
    quantifier_eliminate,
)


def test_universal_qe_fiber_sector_sample_is_interior():
    t, e = sp.symbols("t e", real=True)
    matrix = sp.Implies(
        sp.And(t > 0, t < e),
        sp.Eq(t**4 + t**2 - t, 0),
    )
    result = quantifier_eliminate(
        matrix,
        (("forall", t),),
        variables=(e, t),
        strategy="complete",
    )
    assert result == (e <= 0)


@pytest.mark.parametrize(
    ("symbol", "predicate"),
    [
        (sp.Symbol("p", positive=True), lambda value: value > 0),
        (sp.Symbol("n", negative=True), lambda value: value < 0),
        (sp.Symbol("z", nonzero=True, real=True), lambda value: value != 0),
    ],
)
def test_find_instance_preserves_symbol_sign_assumptions(symbol, predicate):
    point = find_instance(sp.true, (symbol,))
    assert point is not None
    assert predicate(point[symbol])


def test_local_components_do_not_merge_branches_that_reconnect_far_away():
    x, y = sp.symbols("x y", real=True)
    branches = sp.Or(
        sp.And(sp.Eq(y, x), x > 0, x <= 1),
        sp.And(sp.Eq(y, -x), x > 0, x <= 1),
        sp.And(sp.Eq(x, 1), y >= -1, y <= 1),
    )
    components = local_components(branches, (0, 0), (x, y))
    assert len(components) == 2


def test_local_preimage_rejects_nonunique_exact_fiber():
    x, y = sp.symbols("x y", real=True)
    target = local_germ(y >= 0, (1,), (y,))
    with pytest.raises(ValueError, match="unique"):
        local_preimage(target, (x**2,), (x,))


def test_parameter_strata_support_custom_invariant_cases():
    x, a = sp.symbols("x a", real=True)
    strata = parameter_strata(
        sp.Eq(x**2, a),
        (x,),
        (a,),
        conditions_by_value={
            "zero-root": sp.Eq(x, 0),
            "nonzero-root": sp.Ne(x, 0),
        },
    )
    assert tuple(s.value for s in strata) == ("zero-root", "nonzero-root")
    assert sp.simplify(strata[0].condition.subs(a, 0)) is sp.true
    assert sp.simplify(strata[1].condition.subs(a, 1)) is sp.true


def test_sign_strata_prune_impossible_prefixes():
    x = sp.Symbol("x", real=True)
    strata = local_sign_strata((x, x**2), sp.true, (0,), (x,))
    assert {s.signs for s in strata} == {(-1, 1), (0, 0), (1, 1)}
