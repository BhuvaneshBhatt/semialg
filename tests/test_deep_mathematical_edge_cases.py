from __future__ import annotations

import pytest
import sympy as sp

from semialg import connected_components, is_connected
from semialg.algebraic.rational_univariate import (
    compute_rational_univariate_representation,
    solve_zero_dimensional_system_with_rur,
)
from semialg.cad_algorithms.reduced import decomp_form_reduced_safe
from semialg.formula import parse_formula
from semialg.parameters import root_count_conditions
from semialg.tticad.safe import decompose_tticad_safe


def _selected_count(conditions, assignment):
    selected = []
    for count, condition in conditions.items():
        value = sp.simplify(condition.subs(assignment))
        if value is sp.true or value == sp.true:
            selected.append(count)
    assert len(selected) == 1
    return selected[0]


def _distinct_real_root_count(expr, variable):
    roots = sp.real_roots(sp.Poly(expr, variable))
    return len(tuple(dict.fromkeys(roots)))


def test_cubic_parameter_classification_handles_degree_drop_exactly():
    x, a, b, c = sp.symbols("x a b c", real=True)
    conditions = root_count_conditions(a * x**3 + b * x + c, x, [a, b, c])

    assert _selected_count(conditions, {a: 0, b: 2, c: 3}) == 1
    assert _selected_count(conditions, {a: 0, b: 0, c: 3}) == 0
    assert _selected_count(conditions, {a: 0, b: 0, c: 0}) is sp.oo


def test_cubic_parameter_strata_match_specialized_exact_roots():
    x, a, b, c, d = sp.symbols("x a b c d", real=True)
    family = a * x**3 + b * x**2 + c * x + d
    conditions = root_count_conditions(family, x, [a, b, c, d])
    cases = (
        (-2, -1, 0, 1),
        (-1, 0, 1, 0),
        (0, 2, -1, 3),
        (0, 0, 2, -4),
        (0, 0, 0, 5),
        (1, -3, 3, -1),
        (1, 0, -1, 0),
        (2, 1, -2, -1),
    )
    for coefficients in cases:
        assignment = dict(zip((a, b, c, d), coefficients, strict=True))
        specialized = sp.expand(family.subs(assignment))
        selected = _selected_count(conditions, assignment)
        if specialized == 0:
            assert selected is sp.oo
        else:
            assert selected == _distinct_real_root_count(specialized, x)


def test_degree_five_parameter_family_is_exactly_sturm_stratified():
    x, a = sp.symbols("x a", real=True)
    conditions = root_count_conditions(x**5 + a * x + 1, x, [a])

    assert sp.Integer(-1) not in conditions
    assert set(conditions) == {sp.Integer(1), sp.Integer(2), sp.Integer(3)}


def test_mccallum_nullification_falls_back_but_lazard_remains_certified():
    x, y = sp.symbols("x y", real=True)
    formula = parse_formula(sp.And(sp.Eq(x * y, 0), y > 0))

    mccallum = decomp_form_reduced_safe(formula, (x, y), backend="mccallum")
    lazard = decomp_form_reduced_safe(formula, (x, y), backend="lazard")

    assert mccallum.complete
    assert mccallum.used_fallback
    assert mccallum.side_conditions is not None
    assert mccallum.side_conditions.nullification_events

    assert lazard.complete
    assert not lazard.used_fallback
    assert lazard.certificate is not None and lazard.certificate.valid
    assert lazard.side_conditions is not None and lazard.side_conditions.valid


def test_tticad_uses_family_aware_safe_driver():
    x, y = sp.symbols("x y", real=True)
    formula = parse_formula(
        sp.Or(
            sp.And(sp.Eq(y**2 - x, 0), y >= 0),
            sp.And(sp.Eq(y**2 + x, 0), y <= 0),
        )
    )

    result = decompose_tticad_safe(formula, (x, y))

    assert result.complete
    assert result.certificate is not None
    if not result.used_fallback:
        assert result.certificate.valid
        assert result.validity.valid
        assert result.projection_validity.valid


@pytest.mark.parametrize(
    "scalar",
    [sp.Integer(1), sp.sqrt(2), sp.sqrt(3), 1 + sp.sqrt(2), sp.sqrt(2) + sp.sqrt(3)],
)
def test_rur_solution_is_invariant_under_nonzero_algebraic_scaling(scalar):
    x, y = sp.symbols("x y", real=True)
    base = (x**2 - 2, y - x)
    equations = tuple(sp.expand(scalar * equation) for equation in base)

    solutions = solve_zero_dimensional_system_with_rur(equations, (x, y), real=True)

    assert solutions == ((-sp.sqrt(2), -sp.sqrt(2)), (sp.sqrt(2), sp.sqrt(2)))


def test_rur_cache_and_domain_transition_do_not_change_solutions():
    x, y = sp.symbols("x y", real=True)
    rational = (x**2 - 2, y - x)
    algebraic = tuple(sp.expand(sp.sqrt(3) * equation) for equation in rational)

    rational_rep = compute_rational_univariate_representation(rational, (x, y))
    algebraic_rep = compute_rational_univariate_representation(algebraic, (x, y))
    rational_again = solve_zero_dimensional_system_with_rur(rational, (x, y), real=True)

    assert str(rational_rep.defining_polynomial.domain) == "QQ"
    assert getattr(algebraic_rep.defining_polynomial.domain, "is_AlgebraicField", False)
    assert rational_again == ((-sp.sqrt(2), -sp.sqrt(2)), (sp.sqrt(2), sp.sqrt(2)))


@pytest.mark.parametrize(
    ("parameter", "component_count"),
    [(-1, 2), (0, 1), (1, 1)],
)
def test_singular_fiber_topology_changes_are_detected(parameter, component_count):
    x, y = sp.symbols("x y", real=True)
    fiber = sp.Eq(y**2, x**2 * (x + parameter))

    components = connected_components(fiber, (x, y))

    assert len(components) == component_count
    assert is_connected(fiber, (x, y)) is (component_count == 1)


def test_reducible_disjoint_circles_are_two_components():
    x, y = sp.symbols("x y", real=True)
    region = sp.Eq((x**2 + y**2 - 1) * (x**2 + y**2 - 4), 0)

    components = connected_components(region, (x, y))

    assert len(components) == 2
    assert not is_connected(region, (x, y))


def test_reducible_intersecting_axes_merge_into_one_component():
    x, y = sp.symbols("x y", real=True)
    region = sp.Eq(x * y, 0)

    assert len(connected_components(region, (x, y))) == 1
    assert is_connected(region, (x, y))
