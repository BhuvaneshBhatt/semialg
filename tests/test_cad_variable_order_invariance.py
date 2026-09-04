import itertools

import pytest
import sympy as sp

from semialg import is_equal
from semialg.formula import parse_formula
from semialg.qe.complete import qe_by_complete_cad

x, y, z, a = sp.symbols("x y z a", real=True)

CASES = [
    (sp.And(x**2 <= a, x >= 0), [("exists", x)], [a]),
    (sp.And(sp.Eq(x**2, a), x >= 0), [("exists", x)], [a]),
    (sp.And(x**2 + x <= a, x >= -1), [("exists", x)], [a]),
    (sp.And(x >= 0, x <= a), [("exists", x)], [a]),
]


def _run(expr, quantifiers, free, variables, strategy):
    return qe_by_complete_cad(
        variables,
        quantifiers,
        parse_formula(expr),
        free_variables=free,
        variable_order_strategy=strategy,
        return_result=True,
    )


@pytest.mark.parametrize("expr,quantifiers,free", CASES)
def test_auto_order_semantically_matches_preserved_order(expr, quantifiers, free):
    qvars = [v for _, v in quantifiers]
    variables = [*free, *qvars]
    auto = _run(expr, quantifiers, free, variables, "auto")
    preserved = _run(expr, quantifiers, free, variables, "preserve")
    assert is_equal(auto.formula, preserved.formula, free)


TWO_VAR_TRUTH_CASES = [
    sp.And(x**2 <= 1, y**2 <= 1, x + y > 0),
    sp.And(sp.Eq(x + y, 1), x >= 0, y >= 0),
    sp.And(x * y == 1, x > 0, y > 0),
    sp.And(x**2 + y**2 < 1, x > 2),
    sp.And(x**2 <= y, y <= 1),
    sp.And((x - y) ** 2 <= 4, x + y == 0),
]


@pytest.mark.parametrize("expr", TWO_VAR_TRUTH_CASES)
@pytest.mark.parametrize("order", [(x, y), (y, x)])
def test_existential_two_variable_orders_give_same_truth(expr, order):
    quantifiers = [("exists", v) for v in order]
    result = _run(expr, quantifiers, [], list(order), "preserve")
    reference = _run(expr, [("exists", x), ("exists", y)], [], [x, y], "preserve")
    assert result.truth_value == reference.truth_value


@pytest.mark.parametrize("strategy", ["auto", "brown", "projection"])
def test_ordering_never_crosses_two_block_quantifier_boundary(strategy):
    expr = x**2 <= y
    result = _run(expr, [("exists", x), ("forall", y)], [], [x, y], strategy)
    assert result.quantifiers == (("exists", x), ("forall", y))


@pytest.mark.parametrize(
    "expr",
    [
        sp.And(x**2 + y**2 <= a, x >= 0),
        sp.And(sp.Eq(x + y, a), x >= 0, y >= 0),
    ],
)
def test_auto_and_projection_scoring_have_identical_qe_semantics(expr):
    quantifiers = [("exists", x), ("exists", y)]
    auto = _run(expr, quantifiers, [a], [a, x, y], "auto")
    projection = _run(expr, quantifiers, [a], [a, x, y], "projection")
    assert is_equal(auto.formula, projection.formula, [a])


@pytest.mark.slow
def test_three_variable_same_block_permutation_invariance():
    expr = sp.And(x**2 <= 1, y**2 <= 1, z**2 <= 1, x + y + z > 0)
    values = []
    for order in itertools.permutations((x, y, z)):
        result = _run(expr, [("exists", v) for v in order], [], list(order), "preserve")
        values.append(result.truth_value)
    assert all(value == values[0] for value in values)


def test_suggest_variable_order_public_api_returns_all_variables_once():
    import sympy as sp

    from semialg.formula import parse_formula
    from semialg.heuristics import suggest_variable_order

    x, y = sp.symbols("x y", real=True)
    order = suggest_variable_order(parse_formula((x**2 + y**2 <= 1) & (x + y >= 0)))

    assert set(order) == {x, y}
    assert len(order) == 2


def test_preserved_same_kind_quantifier_permutations_share_canonical_cad():
    from semialg.cad_algorithms import cad_cache_stats, clear_cad_caches

    expr = sp.And(x**2 <= 1, y**2 <= 1, x + y > 0)
    clear_cad_caches()
    first = _run(expr, [("exists", y), ("exists", x)], [], [y, x], "preserve")
    after_first = cad_cache_stats()
    second = _run(expr, [("exists", x), ("exists", y)], [], [x, y], "preserve")
    after_second = cad_cache_stats()

    assert first.truth_value == second.truth_value
    assert first.quantifiers == (("exists", x), ("exists", y))
    assert second.quantifiers == first.quantifiers
    assert after_second.complete_cad_hits > after_first.complete_cad_hits
