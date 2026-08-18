import sympy as sp

from semialg.algebraic.rational_univariate import (
    compute_rational_univariate_representation,
    solve_formula_with_rur,
    solve_rur_representation,
    solve_zero_dimensional_system_with_rur,
)
from semialg.qe.virtual_substitution.witness import _ordered_real_values
from semialg.solve.reduce import reduce_text


def test_partial_rur_disjunction_is_unknown_not_unsat():
    x = sp.symbols("x")
    result = solve_formula_with_rur(sp.Or(sp.Eq(x**2 + 1, 0), x > 0), [x], real=True)
    assert result is not None
    assert result.status == "unknown"
    assert result.skipped_branches == 1
    assert not result.complete


def test_reduce_does_not_use_partial_rur_as_false():
    x = sp.symbols("x")
    result = reduce_text(
        "exists x. (x^2 + 1 == 0) | (x > 0)",
        variable_order=[x],
        strategy="auto",
        return_result=True,
    )
    assert result.result == sp.true
    assert result.method != "rational_univariate"


def test_solve_rur_representation_reuses_existing_representation():
    x, y = sp.symbols("x y")
    rep = compute_rational_univariate_representation([x**2 + y**2 - 1, x - y], [x, y])
    assert solve_rur_representation(rep, real=True) == solve_zero_dimensional_system_with_rur(
        [x**2 + y**2 - 1, x - y], [x, y]
    )


def test_complex_rur_uses_crootof_for_irreducible_cubic():
    x = sp.symbols("x")
    roots = solve_zero_dimensional_system_with_rur([x**3 - x + 1], [x], real=False)
    assert len(roots) == 3
    assert all(isinstance(root[0], sp.CRootOf) for root in roots)


def test_virtual_substitution_orders_algebraic_roots_exactly():
    ordered = _ordered_real_values([sp.sqrt(2), sp.Rational(1, 3), -sp.sqrt(2)])
    assert ordered == [-sp.sqrt(2), sp.Rational(1, 3), sp.sqrt(2)]
