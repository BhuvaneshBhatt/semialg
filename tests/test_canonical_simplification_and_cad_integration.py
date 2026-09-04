import sympy as sp
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg import integrate_over_region, reduce_region_integral, simplify_boole, simplify_region
from semialg.simplify import simplify_semialgebraic_formula

x, y, z = sp.symbols("x y z", real=True)


def test_polynomial_relations_have_stable_scalar_and_orientation_normal_form():
    forms = (
        2 * x - 2 * y > 0,
        y - x < 0,
        7 * x - 7 * y > 0,
    )
    simplified = tuple(
        simplify_semialgebraic_formula(form, implication_minimize=False) for form in forms
    )
    assert simplified == (x - y > 0,) * 3


def test_zero_set_multiplicity_is_removed_canonically():
    assert simplify_semialgebraic_formula(
        sp.Eq(6 * (x - 1) ** 4 * (x + 2) ** 2, 0), implication_minimize=False
    ) == sp.Eq(x**2 + x - 2, 0)
    assert simplify_semialgebraic_formula(
        sp.Ne(-3 * (x - 1) ** 6, 0), implication_minimize=False
    ) == sp.Ne(x, 1)


def test_canonical_simplification_is_idempotent_and_shared_by_public_simplifiers():
    expr = sp.Or(
        sp.And(2 * x - 2 * y > 0, x > -1),
        sp.And(y - x < 0, x > -1),
    )
    first = simplify_semialgebraic_formula(expr)
    second = simplify_semialgebraic_formula(first)
    assert first == second == sp.And(x > -1, x - y > 0)
    assert simplify_boole(expr, [x, y]) == first
    assert simplify_region(expr, [x, y]).formula == first


@given(st.integers(min_value=1, max_value=20))
@settings(max_examples=20, deadline=None)
def test_positive_integer_scaling_does_not_change_canonical_inequality(scale):
    expected = simplify_semialgebraic_formula(x**2 + y - 1 <= 0, implication_minimize=False)
    actual = simplify_semialgebraic_formula(scale * (x**2 + y - 1) <= 0, implication_minimize=False)
    assert actual == expected


def test_cad_integration_selects_polynomial_coordinate_order_in_2d():
    condition = sp.And(y >= 0, y <= 1, x >= y**2, x <= y)
    reduced = reduce_region_integral(1, condition, [x, y])

    assert reduced.method == "cad_variable_order_cell_integration"
    assert len(reduced.pieces) == 1
    piece = reduced.pieces[0]
    assert piece.limits == ((x, y**2, y), (y, 0, 1))
    assert piece.diagnostics["integration_variable_order"] == (y, x)
    assert piece.diagnostics["typed_bounds_verified"] is True
    assert integrate_over_region(1, condition, [x, y]) == sp.Rational(1, 6)


def test_cad_integration_order_search_extends_to_nested_3d_regions():
    condition = sp.And(
        z >= 0,
        z <= 1,
        y >= z**2,
        y <= z,
        x >= y**2,
        x <= y,
    )
    reduced = reduce_region_integral(1, condition, [x, y, z])

    assert reduced.method == "cad_variable_order_cell_integration"
    assert len(reduced.pieces) == 1
    piece = reduced.pieces[0]
    assert piece.limits == ((x, y**2, y), (y, z**2, z), (z, 0, 1))
    assert piece.diagnostics["integration_variable_order"] == (z, y, x)
    assert piece.diagnostics["typed_bounds_verified"] is True
    assert integrate_over_region(1, condition, [x, y, z]) == sp.Rational(13, 420)
