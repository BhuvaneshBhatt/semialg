import pytest
import sympy as sp

import semialg


def test_equivalent_does_not_conflate_same_name_symbols_with_different_assumptions():
    x = sp.Symbol("x")
    x_real = sp.Symbol("x", real=True)

    assert semialg.equivalent(x > 0, x_real > 0) is False


def test_equivalent_accepts_generator_variable_sequence_without_consuming_it():
    x = sp.Symbol("x", real=True)
    variables = (var for var in (x,))

    assert semialg.equivalent(x > 0, x >= 0, variables) is False


def test_measure_rejects_bound_for_undeclared_symbol():
    x, y = sp.symbols("x y", real=True)

    with pytest.raises(ValueError, match="bound variable"):
        semialg.semialgebraic_measure(x > 0, [x], bounds={y: (0, 1)})


def test_integration_rejects_unknown_string_bound():
    x = sp.Symbol("x", real=True)

    with pytest.raises(ValueError, match="bound variable"):
        semialg.integrate_over_region(1, x > 0, [x], bounds={"y": (0, 1)})


def test_measure_rejects_reversed_bounds_instead_of_returning_negative_measure():
    x = sp.Symbol("x", real=True)

    with pytest.raises(ValueError, match="lower bound exceeds"):
        semialg.semialgebraic_measure(sp.true, [x], bounds={x: (2, 1)})


def test_bounds_mapping_requires_lower_upper_pair():
    x = sp.Symbol("x", real=True)

    with pytest.raises(ValueError, match="lower, upper"):
        semialg.semialgebraic_measure(sp.true, [x], bounds={x: (0,)})


def test_bounds_sequence_rejects_duplicate_variable():
    x = sp.Symbol("x", real=True)

    with pytest.raises(ValueError, match="duplicate bound"):
        semialg.semialgebraic_measure(
            sp.true,
            [x],
            bounds=[(x, 0, 1), (x, 0, 2)],
        )


def test_standard_region_integration_preserves_input_symbol_identity():
    x = sp.Symbol("x")

    result = semialg.integrate_over_standard_region(
        x,
        semialg.IntervalRegion(0, 1),
        ["x"],
    )

    assert result == sp.Rational(1, 2)


def test_standard_region_intersection_orders_close_endpoints_exactly():
    x = sp.Symbol("x", real=True)
    epsilon = sp.Rational(1, 10**120)
    region = semialg.BooleanRegion(
        "intersection",
        (
            semialg.IntervalRegion(0, 1 + epsilon),
            semialg.IntervalRegion(0, 1),
        ),
    )

    assert semialg.integrate_over_standard_region(1, region, [x]) == 1


def test_region_components_does_not_shrink_contained_interval():
    x = sp.Symbol("x", real=True)
    region = sp.Or(
        sp.And(x >= 0, x <= 10),
        sp.And(x >= 1, x <= 2),
        evaluate=False,
    )

    components = semialg.region_components(region, [x])

    assert len(components) == 1
    assert semialg.equivalent(components[0], sp.And(x >= 0, x <= 10), [x])


def test_region_components_preserves_strictness_at_shared_endpoint():
    x = sp.Symbol("x", real=True)
    region = sp.Or(sp.And(x > 0, x < 1), sp.And(x > 1, x < 2))

    components = semialg.region_components(region, [x])

    assert len(components) == 2
