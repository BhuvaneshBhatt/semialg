from __future__ import annotations

import sympy as sp

import semialg
from semialg.affine_reduction import parametric_affine_reduction
from semialg.context import SemialgebraicContext
from semialg.matrix_analysis import (
    constant_symmetric_inertia,
    matrix_pd_on,
    matrix_psd_on,
    matrix_rank_on,
    matrix_rank_stratification,
)
from semialg.strict_feasibility import affine_relative_interior_formula, strict_feasible


def test_constant_matrix_inertia_and_definiteness():
    matrix = sp.Matrix([[2, 1], [1, 2]])
    assert constant_symmetric_inertia(matrix) == (2, 0, 0)
    assert matrix_psd_on(matrix)
    assert matrix_pd_on(matrix)
    assert semialg.matrix_definiteness(matrix, requested="positive_definite")


def test_polynomial_matrix_psd_returns_counterexample_point():
    x = sp.symbols("x", real=True)
    result = matrix_psd_on(
        sp.Matrix([[x]]), [x], domain=sp.And(x >= -1, x <= 1), return_result=True
    )
    assert result.outcome is False
    assert result.counterexample is not None
    assert sp.simplify(result.counterexample[x]) < 0
    assert result.counterexample_vector == (1,)


def test_parameterized_matrix_psd_condition():
    a = sp.symbols("a", real=True)
    result = matrix_psd_on(sp.Matrix([[a]]), [], parameters=[a])
    assert result.select({a: 2}) is True
    assert result.select({a: -2}) is False


def test_matrix_rank_stratification_and_constant_rank_on_region():
    a, x = sp.symbols("a x", real=True)
    strata = matrix_rank_stratification(sp.Matrix([[a, 0], [0, 1]]), [a])
    assert strata.select({a: 0}) == 1
    assert strata.select({a: 3}) == 2
    assert matrix_rank_on(sp.Matrix([[1, x], [0, 0]]), [x]) == 1


def test_relative_strict_feasibility_preserves_affine_hull_equalities():
    x, y = sp.symbols("x y", real=True)
    constraints = sp.And(sp.Eq(y, 0), x >= 0, x <= 1)
    relative_formula = affine_relative_interior_formula(constraints, [x, y])
    assert relative_formula.has(sp.Eq(y, 0))
    result = strict_feasible(constraints, [x, y], return_result=True)
    assert result.feasible
    assert result.witness is not None
    assert result.witness[y] == 0
    assert 0 < result.witness[x] < 1


def test_relative_feasibility_preserves_hull_equality():
    x, y = sp.symbols("x y", real=True)
    constraints = sp.And(sp.Eq(y, 0), y >= 0, x >= 0, x <= 1)
    assert strict_feasible(constraints, [x, y], relative=True)


def test_sign_primitives_and_parameter_conditions():
    x, a = sp.symbols("x a", real=True)
    assert semialg.prove_zero(x - x, [x])
    assert semialg.prove_nonzero(x**2 + 1, [x])
    assert semialg.function_sign(x**2, [x]) == "nonnegative"
    parameter_result = semialg.prove_nonnegative(
        x + a, [x], assumptions=sp.And(x >= 0, x <= 1), parameters=[a]
    )
    assert parameter_result.select({a: 0}) is True
    assert parameter_result.select({a: -1}) is False


def test_parametric_affine_reduction_branches_on_pivot_exception():
    x, a = sp.symbols("x a", real=True)
    result = parametric_affine_reduction(sp.Eq(a * x + 1, 0), [x], [a])
    generic = result.select({a: 2})
    exceptional = result.select({a: 0})
    assert generic.substitution_map[x] == -1 / a
    assert generic.formula is sp.true or generic.formula == sp.true
    assert exceptional.formula is sp.false or exceptional.formula == sp.false


def test_semialgebraic_context_reuses_new_reasoning_interfaces():
    x = sp.symbols("x", real=True)
    context = SemialgebraicContext(sp.And(x >= 0, x <= 1), (x,))
    assert context.function_sign(x) == "nonnegative"
    assert context.matrix_definiteness(sp.Matrix([[x]]))
    assert context.strict_feasible()
    derived = context.add_constraints(x <= sp.Rational(1, 2))
    assert derived.computation is context.computation
    assert derived.implies(x <= 1)


def test_polynomial_convexity_reports_matrix_backend():
    import semialg.convexity as convexity_module

    x = sp.symbols("x", real=True)
    certificate = convexity_module.polynomial_convexity_certificate(x**4 + x**2, [x])
    assert certificate.outcome is True
    assert certificate.method == "principal-minor-shared-cad"


def test_function_sign_intersects_natural_real_domain():
    x = sp.symbols("x", real=True)
    assert semialg.function_sign(sp.sqrt(x), [x], assumptions=x < 0) == "empty_domain"
    assert semialg.function_sign(sp.sqrt(x), [x]) == "nonnegative"
