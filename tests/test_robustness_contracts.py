from __future__ import annotations

import sympy as sp

from semialg import (
    analyze_affine_map,
    function_mapping_properties,
    is_bijective,
    is_function_continuous,
    is_function_smooth,
    is_injective,
    is_surjective,
    matrix_definiteness,
    polynomial_nonnegative_decision,
    real_algebraic_feasibility,
    tangent_cone,
    tangent_space,
)
from semialg.decision_portfolio import CertifiedDecisionResult
from semialg.real_algebraic import RealAlgebraicFeasibilityResult


def test_explicit_tangent_variables_do_not_absorb_parameter_named_d0():
    x = sp.Symbol("x", real=True)
    d0 = sp.Symbol("d0", real=True)

    space = tangent_space(x - d0, {x: d0}, (x,))
    cone = tangent_cone(x - d0, {x: d0}, (x,))

    assert space.variables == (x,)
    assert cone.variables == (x,)
    assert all(d0 not in equation.free_symbols for equation in space.equations)
    assert all(direction.name != "d0" for direction in cone.direction_variables)


def test_image_coordinates_avoid_parameter_named_y0():
    x = sp.Symbol("x", real=True)
    y0 = sp.Symbol("y0", real=True)
    props = function_mapping_properties(x + y0, (x,))
    assert props.image_variables[0].name != "y0"


def test_parameterized_affine_determinant_is_not_claimed_everywhere_nonzero():
    x = sp.Symbol("x", real=True)
    a = sp.Symbol("a", real=True)
    analysis = analyze_affine_map((a * x,), (x,))
    assert analysis.determinant == a
    assert analysis.invertible is None


def test_constant_matrix_unknown_remains_unknown():
    unknown_constant = sp.Function("c")(1)
    assert matrix_definiteness([[unknown_constant]], requested="positive_definite") is None
    detailed = matrix_definiteness(
        [[unknown_constant]], requested="positive_definite", return_result=True
    )
    assert detailed.outcome is None
    assert detailed.certified is False


def test_single_answer_decisions_default_to_mathematical_values():
    x = sp.Symbol("x", real=True)
    assert polynomial_nonnegative_decision(x**2 + 1, (x,), sos_backend="none") is True
    detailed = polynomial_nonnegative_decision(
        x**2 + 1, (x,), sos_backend="none", return_result=True
    )
    assert isinstance(detailed, CertifiedDecisionResult)
    assert detailed.decision is True

    assert real_algebraic_feasibility((x**2 - 1,), (x,)) is True
    ars = real_algebraic_feasibility((x**2 - 1,), (x,), return_result=True)
    assert isinstance(ars, RealAlgebraicFeasibilityResult)
    assert ars.satisfiable is True


def test_function_property_convenience_predicates_project_aggregate_results():
    x = sp.Symbol("x", real=True)
    assert is_injective(x, x) is True
    assert is_surjective(x, x) is True
    assert is_bijective(x, x) is True
    assert is_function_continuous(x**2, x) is True
    assert is_function_smooth(x**2, x) is True
