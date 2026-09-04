from __future__ import annotations

import pytest
import sympy as sp
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg import cad, is_equal
from semialg.formula import parse_formula
from semialg.qe.complete import qe_by_complete_cad


@st.composite
def _finite_triangular_cases(draw):
    radius = draw(st.integers(min_value=1, max_value=3))
    slope = draw(st.integers(min_value=-2, max_value=2).filter(lambda value: value != 0))
    shift = draw(st.integers(min_value=-2, max_value=2))
    threshold = draw(st.integers(min_value=-4, max_value=4))
    return radius, slope, shift, threshold


def _formula(case):
    radius, slope, shift, threshold = case
    x, y = sp.symbols("x y", real=True)
    return (x, y), sp.And(
        sp.Eq(x**2 - radius**2, 0),
        sp.Eq(y - slope * x - shift, 0),
        y > threshold,
    )


@settings(max_examples=16, deadline=None)
@given(_finite_triangular_cases())
def test_auto_variety_cad_matches_forced_collins(case) -> None:
    variables, formula = _formula(case)
    auto = cad(formula, variables)
    collins = cad(formula, variables, strategy="collins")
    assert is_equal(auto, collins, variables)


@pytest.mark.slow
@settings(max_examples=80, deadline=None)
@given(_finite_triangular_cases())
def test_large_auto_variety_cad_matches_forced_collins(case) -> None:
    variables, formula = _formula(case)
    auto = cad(formula, variables)
    collins = cad(formula, variables, strategy="collins")
    assert is_equal(auto, collins, variables)


@settings(max_examples=10, deadline=None)
@given(_finite_triangular_cases())
def test_existential_qe_variety_matches_forced_collins(case) -> None:
    (x, y), formula = _formula(case)
    matrix = parse_formula(formula)
    common = dict(
        vars_=(x, y),
        quantifiers=(("exists", y),),
        matrix=matrix,
        free_variables=(x,),
        variable_order_strategy="preserve",
        use_presolve=False,
        return_result=True,
    )
    variety = qe_by_complete_cad(allow_variety_cad=True, **common)
    collins = qe_by_complete_cad(allow_variety_cad=False, **common)
    assert is_equal(variety.formula, collins.formula, (x,))


@pytest.mark.slow
@settings(max_examples=120, deadline=None)
@given(_finite_triangular_cases())
def test_large_existential_qe_variety_matches_forced_collins(case) -> None:
    (x, y), formula = _formula(case)
    matrix = parse_formula(formula)
    common = dict(
        vars_=(x, y),
        quantifiers=(("exists", y),),
        matrix=matrix,
        free_variables=(x,),
        variable_order_strategy="preserve",
        use_presolve=False,
        return_result=True,
    )
    variety = qe_by_complete_cad(allow_variety_cad=True, **common)
    collins = qe_by_complete_cad(allow_variety_cad=False, **common)
    assert is_equal(variety.formula, collins.formula, (x,))


@st.composite
def _compatible_root_cases(draw):
    x_radius = draw(st.integers(min_value=1, max_value=3))
    y_radius = draw(st.integers(min_value=1, max_value=3))
    orientation = draw(st.sampled_from((-1, 1)))
    threshold = draw(st.integers(min_value=-4, max_value=4))
    return x_radius, y_radius, orientation, threshold


def _compatible_formula(case):
    x_radius, y_radius, orientation, threshold = case
    x, y = sp.symbols("x y", real=True)
    relation = x * y - orientation * x_radius * y_radius
    return (x, y), sp.And(
        sp.Eq(x**2 - x_radius**2, 0),
        sp.Eq(y**2 - y_radius**2, 0),
        sp.Eq(relation, 0),
        x + y > threshold,
    )


@settings(max_examples=16, deadline=None)
@given(_compatible_root_cases())
def test_compatible_root_variety_matches_forced_collins(case) -> None:
    variables, formula = _compatible_formula(case)
    auto = cad(formula, variables)
    collins = cad(formula, variables, strategy="collins")
    assert is_equal(auto, collins, variables)


@pytest.mark.slow
@settings(max_examples=24, deadline=None)
@given(_compatible_root_cases())
def test_large_compatible_root_variety_matches_forced_collins(case) -> None:
    variables, formula = _compatible_formula(case)
    auto = cad(formula, variables)
    collins = cad(formula, variables, strategy="collins")
    assert is_equal(auto, collins, variables)


@st.composite
def _scaled_generator_cases(draw):
    root = draw(st.integers(min_value=1, max_value=3))
    slope = draw(st.sampled_from((-2, -1, 1, 2)))
    shift = draw(st.integers(min_value=-2, max_value=2))
    scales = draw(
        st.tuples(
            st.sampled_from((-3, -2, -1, 1, 2, 3)),
            st.sampled_from((-3, -2, -1, 1, 2, 3)),
        )
    )
    reverse = draw(st.booleans())
    swap_order = draw(st.booleans())
    return root, slope, shift, scales, reverse, swap_order


@settings(max_examples=16, deadline=None)
@given(_scaled_generator_cases())
def test_generator_scaling_and_variable_order_are_semantic_invariants(case) -> None:
    root, slope, shift, scales, reverse, swap_order = case
    x, y = sp.symbols("x y", real=True)
    generators = [x**2 - root**2, y - slope * x - shift]
    reference = cad(sp.And(*(sp.Eq(g, 0) for g in generators)), (x, y))

    transformed = [scale * generator for scale, generator in zip(scales, generators, strict=True)]
    if reverse:
        transformed.reverse()
    order = (y, x) if swap_order else (x, y)
    candidate = cad(sp.And(*(sp.Eq(g, 0) for g in transformed)), order)

    assert is_equal(reference, candidate, (x, y))


@st.composite
def _common_boolean_cases(draw):
    root = draw(st.integers(min_value=1, max_value=3))
    scale = draw(st.sampled_from((-4, -3, -2, -1, 1, 2, 3, 4)))
    threshold = draw(st.integers(min_value=-3, max_value=3))
    return root, scale, threshold


@settings(max_examples=12, deadline=None)
@given(_common_boolean_cases())
def test_nested_boolean_common_equalities_match_collins(case) -> None:
    root, scale, threshold = case
    x, y = sp.symbols("x y", real=True)
    common = x**2 - root**2
    formula = sp.And(
        sp.Or(
            sp.And(sp.Eq(common, 0), sp.Eq(y - x, 0), y > threshold),
            sp.And(sp.Eq(scale * common, 0), sp.Eq(y + x, 0), y < threshold),
        ),
        x != 0,
    )

    auto = cad(formula, (x, y))
    collins = cad(formula, (x, y), strategy="collins")
    assert is_equal(auto, collins, (x, y))
