"""Reusable differential checks compare mathematical semantics, not representation shape."""

from __future__ import annotations

import sympy as sp

from semialg import is_tautology, polynomial_nonnegative_decision
from semialg.algebraic.rational_univariate import solve_zero_dimensional_system_with_rur
from semialg.formula import parse_formula
from semialg.preprocess import semialgebraicize
from semialg.qe.complete import qe_by_complete_cad

from ._semantic_equivalence import assert_exact_point_sets_equal, assert_formulas_equivalent


def test_groebner_variety_qe_matches_forced_collins_qe_semantically():
    x, y = sp.symbols("x y", real=True)
    matrix = parse_formula(sp.And(sp.Eq(x**2 - 2, 0), sp.Eq(y - x, 0), y > 0))

    variety = qe_by_complete_cad(
        (x, y),
        (("exists", y),),
        matrix,
        free_variables=(x,),
        variable_order_strategy="preserve",
        use_presolve=False,
        allow_variety_cad=True,
        return_result=True,
    )
    collins = qe_by_complete_cad(
        (x, y),
        (("exists", y),),
        matrix,
        free_variables=(x,),
        variable_order_strategy="preserve",
        use_presolve=False,
        allow_variety_cad=False,
        return_result=True,
    )

    assert_formulas_equivalent(variety.formula, collins.formula, (x,))


def test_rur_solution_set_matches_direct_exact_solution_set():
    x, y = sp.symbols("x y", real=True)
    equations = (x**2 - 2, y - x)
    rur = solve_zero_dimensional_system_with_rur(equations, (x, y), real=True)
    direct = tuple((root, root) for root in sp.real_roots(x**2 - 2))

    assert_exact_point_sets_equal(rur, direct, (x, y))


def test_specialized_nonnegativity_decision_matches_complete_formula_truth():
    x, y = sp.symbols("x y", real=True)
    polynomial = x**4 + y**4 + x**2 * y**2 + 1

    specialized = polynomial_nonnegative_decision(polynomial, (x, y), sos_backend="none")
    complete = is_tautology(polynomial >= 0, (x, y))

    assert specialized is complete is True


def test_exact_algebraization_matches_explicit_polynomial_region():
    x = sp.Symbol("x", real=True)
    preprocessed = semialgebraicize(sp.Abs(x) <= 1, variables=(x,))
    (auxiliary,) = preprocessed.aux_vars

    eliminated = qe_by_complete_cad(
        (x, auxiliary),
        (("exists", auxiliary),),
        parse_formula(preprocessed.sympy_expr),
        free_variables=(x,),
        variable_order_strategy="preserve",
        use_presolve=False,
        return_result=True,
    )

    assert_formulas_equivalent(eliminated.formula, x**2 <= 1, (x,))
