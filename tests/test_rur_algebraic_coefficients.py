"""Exact algebraic-coefficient and fallback tests for the RUR backend."""

from __future__ import annotations

import sympy as sp
from sympy.polys.polyerrors import CoercionFailed

from semialg import is_satisfiable
from semialg.algebraic.rational_univariate import (
    compute_rational_univariate_representation,
    solve_formula_with_rur,
    solve_zero_dimensional_system_with_rur,
)


def test_univariate_rur_supports_sqrt2_coefficient_field():
    x = sp.Symbol("x", real=True)

    representation = compute_rational_univariate_representation(
        [x**2 - sp.sqrt(2)],
        [x],
    )
    solutions = solve_zero_dimensional_system_with_rur(
        [x**2 - sp.sqrt(2)],
        [x],
    )

    assert str(representation.defining_polynomial.domain) == "QQ<sqrt(2)>"
    assert solutions == ((-(2 ** sp.Rational(1, 4)),), (2 ** sp.Rational(1, 4),))


def test_multivariate_rur_supports_algebraic_coefficients():
    x, y = sp.symbols("x y", real=True)

    representation = compute_rational_univariate_representation(
        [x + sp.sqrt(2) * y, y - sp.sqrt(2)],
        [x, y],
    )
    solutions = solve_zero_dimensional_system_with_rur(
        [x + sp.sqrt(2) * y, y - sp.sqrt(2)],
        [x, y],
    )

    assert getattr(representation.defining_polynomial.domain, "is_AlgebraicField", False)
    assert solutions == ((-2, sp.sqrt(2)),)


def test_rur_constructs_common_field_for_multiple_algebraic_generators():
    x, y = sp.symbols("x y", real=True)

    representation = compute_rational_univariate_representation(
        [x - sp.sqrt(2), y - sp.sqrt(3)],
        [x, y],
    )
    solutions = solve_zero_dimensional_system_with_rur(
        [x - sp.sqrt(2), y - sp.sqrt(3)],
        [x, y],
    )

    assert getattr(representation.defining_polynomial.domain, "is_AlgebraicField", False)
    assert solutions == ((sp.sqrt(2), sp.sqrt(3)),)


def test_formula_solver_uses_rur_for_algebraic_coefficients():
    x = sp.Symbol("x", real=True)
    formula = sp.Eq(x**2, sp.sqrt(2)) & (x > 0)

    result = solve_formula_with_rur(formula, [x])

    assert result is not None
    assert result.status == "satisfied"
    assert result.assignments == ({x: 2 ** sp.Rational(1, 4)},)


def test_satisfiability_reports_rur_backend_for_algebraic_coefficients():
    x = sp.Symbol("x", real=True)
    formula = sp.Eq(x**2, sp.sqrt(2)) & (x > 0)

    result = is_satisfiable(formula, [x], return_result=True)

    assert result.satisfiable is True
    assert result.method == "rational_univariate"
    assert result.witness == {x: 2 ** sp.Rational(1, 4)}


def test_formula_rur_declines_transcendental_coefficients_cleanly():
    x = sp.Symbol("x", real=True)

    assert solve_formula_with_rur(sp.Eq(x, sp.pi), [x]) is None


def test_formula_rur_declines_symbolic_parameter_coefficients_cleanly():
    x, a = sp.symbols("x a", real=True)

    assert solve_formula_with_rur(sp.Eq(x, a), [x]) is None


def test_formula_rur_contains_raw_coercion_failure():
    import semialg.algebraic.rational_univariate.formula as formula_module

    x = sp.Symbol("x", real=True)

    def fail(*args, **kwargs):
        raise CoercionFailed("synthetic RUR domain failure")

    result = formula_module._solve_formula_with_rur(
        sp.Eq(x, 1),
        [x],
        real=True,
        max_solutions=None,
        branch_solver=fail,
    )

    assert result is None
