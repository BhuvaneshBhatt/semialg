from __future__ import annotations

import itertools

import pytest
import sympy as sp

from semialg import cad, is_equal
from semialg.algebraic.equality_ideal import EqualityIdealContext
from semialg.cad_algorithms.decomposition import (
    decomp_groebner_variety,
    try_decomp_groebner_variety,
)
from semialg.cad_algorithms.lifting.groebner import GroebnerLiftOps
from semialg.cad_algorithms.projection.groebner import build_groebner_variety_projection
from semialg.errors import VarietyCADUnsupported
from semialg.exact_arithmetic import exact_truth
from semialg.formula import equational_constraints, parse_formula
from semialg.presolve import presolve_semialgebraic
from semialg.qe.complete import qe_by_complete_cad


def _raise_sign(_poly, _sample):
    raise ValueError("deliberate exact-sign refusal")


def _raise_fiber(_expr, _variable):
    raise sp.PolynomialError("deliberate fiber refusal")


def test_uncertain_sign_declines_variety_backend_instead_of_pruning_root() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = parse_formula(sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y - x, 0)))
    projection = build_groebner_variety_projection(formula, (x, y))
    assert projection is not None

    ops = GroebnerLiftOps(sign_evaluator=_raise_sign)
    with pytest.raises(VarietyCADUnsupported):
        decomp_groebner_variety(projection, lift_ops=ops)
    assert try_decomp_groebner_variety(projection, lift_ops=ops) is None


def test_uncertain_fiber_conversion_declines_variety_backend() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = parse_formula(sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y - x, 0)))
    projection = build_groebner_variety_projection(formula, (x, y))
    assert projection is not None

    ops = GroebnerLiftOps(fiber_builder=_raise_fiber)
    assert try_decomp_groebner_variety(projection, lift_ops=ops) is None


def test_zero_variable_cad_and_presolve_are_exact() -> None:
    result = cad(sp.true, (), return_result=True)
    assert result.formula is sp.true
    assert result.variables == ()
    assert result.cad.cell_count_by_level() == {0: 1}

    presolved = presolve_semialgebraic(
        sp.Eq(sp.sqrt(2), sp.sqrt(2), evaluate=False), (), eliminate=()
    )
    assert presolved.variables == ()
    assert presolved.linear is True


def test_algebraic_collins_full_elimination() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x**2 - 2, 0), sp.Eq(y - x, 0), y > 0)

    auto = cad(formula, (x, y))
    collins = cad(formula, (x, y), strategy="collins")

    assert is_equal(auto, collins, (x, y))


def test_common_equations_are_canonical_under_scaling_and_sign() -> None:
    x, y = sp.symbols("x y", real=True)
    parsed = parse_formula(
        sp.Or(
            sp.And(sp.Eq(x + y, 0), x > 0),
            sp.And(sp.Eq(-2 * x - 2 * y, 0), y < 0),
        )
    )
    assert equational_constraints(parsed) == [x + y]


def test_scaled_common_equalities_enable_variety_cad_across_boolean_branches() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = sp.Or(
        sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y - x, 0), x > 0),
        sp.And(sp.Eq(2 - 2 * x**2, 0), sp.Eq(3 * y - 3 * x, 0), x < 0),
    )

    result = cad(formula, (x, y), return_result=True)
    assert result.cad.backend == "groebner-variety"
    assert {cell.sample_exprs for cell in result.cad.cells} == {(-1, -1), (1, 1)}


def test_nonradical_quotient_dimension_is_not_geometric_point_count() -> None:
    x, y = sp.symbols("x y", real=True)
    context = EqualityIdealContext((x**2, y - x), (x, y))
    result = cad(sp.And(sp.Eq(x**2, 0), sp.Eq(y - x, 0)), (x, y), return_result=True)

    assert context.dimension == 0
    assert context.quotient_dimension == 2
    assert len(result.cad.cells) == 1
    assert result.cad.cells[0].sample_exprs == (0, 0)


def test_multivariate_complex_only_variety_has_no_real_leaves() -> None:
    x, y = sp.symbols("x y", real=True)
    result = cad(
        sp.And(sp.Eq(x**2 + 1, 0, evaluate=False), sp.Eq(y - x, 0)),
        (x, y),
        return_result=True,
    )
    assert result.cad.backend == "groebner-variety"
    assert result.cad.cells == ()
    assert result.formula is sp.false


def test_algebraic_coefficients_are_supported_by_variety_lifting() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y - sp.sqrt(2) * x, 0), y > 0)
    result = cad(formula, (x, y), return_result=True)

    assert result.cad.backend == "groebner-variety"
    assert len(result.cad.cells) == 1
    point = result.cad.cells[0].sample_exprs
    assert point[0] == 1
    assert sp.simplify(point[1] - sp.sqrt(2)) == 0


def test_variety_leaf_invariants_and_parent_chains() -> None:
    x, y = sp.symbols("x y", real=True)
    equations = (x**2 - 1, y**2 - 1, x * y - 1)
    formula = sp.And(*(sp.Eq(eq, 0) for eq in equations))
    result = cad(formula, (x, y), return_result=True)

    quotient_dimension = result.diagnostics["quotient_dimension"]
    assert quotient_dimension is not None
    assert len(result.cad.cells) <= quotient_dimension
    level_one = {cell.index for cell in result.cad.cells_by_level[1]}
    for leaf in result.cad.cells:
        assert leaf.parent_index in level_one
        substitutions = dict(zip((x, y), leaf.sample_exprs, strict=True))
        assert all(exact_truth(sp.Eq(eq.subs(substitutions), 0)) for eq in equations)


def test_variable_order_permutations_preserve_finite_solution_set() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x**2 - 2, 0), sp.Eq(y - x, 0), y < 0)
    reference = cad(formula, (x, y))
    permuted = cad(formula, (y, x))
    assert is_equal(reference, permuted, (x, y))


def test_generator_permutation_and_nonzero_scaling_preserve_solution_set() -> None:
    x, y = sp.symbols("x y", real=True)
    generators = [x**2 - 1, y - x, x * y - 1]
    reference = cad(sp.And(*(sp.Eq(g, 0) for g in generators)), (x, y))

    for permutation in itertools.permutations(generators):
        scaled = [
            factor * generator for factor, generator in zip((2, -3, 5), permutation, strict=True)
        ]
        candidate = cad(sp.And(*(sp.Eq(g, 0) for g in scaled)), (x, y))
        assert is_equal(reference, candidate, (x, y))


def test_equality_ideal_caches_expensive_queries() -> None:
    x, y = sp.symbols("x y", real=True)
    context = EqualityIdealContext((x**2 - 1, y - x), (x, y))

    initial = context.diagnostics()
    assert initial["standard_materializations"] == 0
    assert context.quotient_dimension == 2
    assert context.diagnostics()["standard_materializations"] == 0

    first_standard = context.standard_exponents
    second_standard = context.standard_exponents
    assert first_standard is second_standard
    assert context.diagnostics()["standard_materializations"] == 1

    assert context.in_radical((x**2 - 1) * (x + 7))
    after_radical = context.diagnostics()["radical_basis_count"]
    assert context.in_radical((x**2 - 1) * (x + 7))
    assert context.diagnostics()["radical_basis_count"] == after_radical

    assert context.has_common_zero_with(x - 1)
    after_common = context.diagnostics()["common_zero_basis_count"]
    assert context.has_common_zero_with(x - 1)
    assert context.diagnostics()["common_zero_basis_count"] == after_common


def test_coordinate_eliminants_are_lazy_in_projection_metadata() -> None:
    x, y, z = sp.symbols("x y z", real=True)
    formula = parse_formula(sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y - x, 0), sp.Eq(z - y, 0)))
    projection = build_groebner_variety_projection(formula, (x, y, z))
    assert projection is not None

    before = projection.equality_context.diagnostics()["fglm_conversion_count"]
    assert before == 1
    assert projection.tower.metadata["coordinate_polynomials"] == "lazy"

    coordinate = projection.coordinate_polynomials
    assert len(coordinate) == 3
    after = projection.equality_context.diagnostics()["fglm_conversion_count"]
    assert after <= 3


def test_existential_qe_can_be_differentially_forced_to_collins() -> None:
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

    assert variety.backend == "groebner-variety-qe"
    assert collins.backend == "collins-complete-qe"
    assert is_equal(variety.formula, collins.formula, (x,))


def test_repeated_and_reducible_finite_variety_keeps_distinct_points_only() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq((x - 1) ** 3 * (x + 1) ** 2, 0), sp.Eq(y - x, 0))
    result = cad(formula, (x, y), return_result=True)

    assert result.cad.backend == "groebner-variety"
    assert {cell.sample_exprs for cell in result.cad.cells} == {(-1, -1), (1, 1)}
    assert result.diagnostics["quotient_dimension"] == 5


def test_rational_coefficients_and_denominators_preserve_variety_dispatch() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(
        sp.Eq(sp.Rational(1, 3) * x**2 - sp.Rational(2, 3), 0),
        sp.Eq(sp.Rational(5, 7) * (y - x), 0),
        y > 0,
    )
    result = cad(formula, (x, y), return_result=True)

    assert result.cad.backend == "groebner-variety"
    assert len(result.cad.cells) == 1
    point = result.cad.cells[0].sample_exprs
    assert sp.simplify(point[0] - sp.sqrt(2)) == 0
    assert sp.simplify(point[1] - sp.sqrt(2)) == 0


def test_boolean_formula_without_common_equalities_uses_complete_space_backend() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = sp.Or(
        sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y - x, 0)),
        sp.And(sp.Eq(x**2 - 4, 0), sp.Eq(y + x, 0)),
    )
    result = cad(formula, (x, y), return_result=True)

    assert result.cad.backend != "groebner-variety"
    assert result.diagnostics["variety_only"] is False


def test_zero_variable_false_formula_has_no_selected_cells() -> None:
    result = cad(sp.false, (), return_result=True)
    assert result.formula is sp.false
    assert result.variables == ()
    assert result.cell_set.cells == ()


def test_common_equation_canonicalization_handles_rational_scaling() -> None:
    x, y = sp.symbols("x y", real=True)
    parsed = parse_formula(
        sp.Or(
            sp.And(sp.Eq(x - 2 * y, 0), x > 0),
            sp.And(sp.Eq(sp.Rational(-3, 5) * x + sp.Rational(6, 5) * y, 0), y < 0),
        )
    )
    assert equational_constraints(parsed) == [x - 2 * y]


def test_quotient_dimension_does_not_materialize_large_standard_basis() -> None:
    variables = sp.symbols("x0:6", real=True)
    context = EqualityIdealContext(tuple(variable**3 for variable in variables), variables)

    assert context.dimension == 0
    assert context.quotient_dimension == 3**6
    assert context.diagnostics()["standard_materializations"] == 0


def test_polynomial_normalization_keeps_algebraic_trailing_variable_factor() -> None:
    from semialg.cad_algorithms.polynomial_utils import normalize_poly

    x, y = sp.symbols("x y", real=True)
    original = sp.Poly(y - sp.sqrt(2) / 2, x, y, extension=True)
    normalized = normalize_poly(original)

    assert normalized is not None
    assert normalized.total_degree() == 1
    assert sp.simplify(normalized.as_expr().subs(y, sp.sqrt(2) / 2)) == 0


def test_standard_monomial_counter_matches_materialized_basis() -> None:
    from semialg.algebraic.rational_univariate.quotient import (
        _standard_exponent_count,
        _standard_exponents,
    )

    cases = (
        (((2, 0), (0, 3)), 2),
        (((3, 0), (1, 1), (0, 2)), 2),
        (((2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0)), 3),
        (((3, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 1)), 3),
    )
    for leading, variable_count in cases:
        materialized = _standard_exponents(leading, variable_count)
        counted = _standard_exponent_count(leading, variable_count)
        assert counted == len(materialized)
