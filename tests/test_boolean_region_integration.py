import sympy as sp

from semialg import integrate_over_region, semialgebraic_measure


def test_measure_disjoint_boolean_union_of_rectangles():
    x, y = sp.symbols("x y", real=True)
    region = sp.Or(
        sp.And(x >= 0, x <= 1, y >= 0, y <= 1),
        sp.And(x >= 2, x <= 4, y >= 0, y <= 1),
    )

    result = semialgebraic_measure(region, [x, y], return_result=True)

    assert sp.simplify(result.value - 3) == 0
    assert result.method == "complete_cad_boolean_cell_integration"
    assert result.diagnostics["boolean_cell_decomposition"] is True


def test_integral_over_overlapping_boolean_union_does_not_double_count_overlap():
    x, y = sp.symbols("x y", real=True)
    region = sp.Or(
        sp.And(x >= 0, x <= 2, y >= 0, y <= 1),
        sp.And(x >= 1, x <= 3, y >= 0, y <= 1),
    )

    value = integrate_over_region(x, region, [x, y])

    assert sp.simplify(value - sp.Rational(9, 2)) == 0


def test_bounded_complement_uses_boolean_cad_cells():
    x, y = sp.symbols("x y", real=True)
    removed = sp.And(x > 0, x < 1, y > 0, y < 1)

    value = semialgebraic_measure(
        sp.Not(removed),
        [x, y],
        bounds={x: (-1, 2), y: (-1, 2)},
    )

    assert sp.simplify(value - 8) == 0


def test_three_dimensional_boolean_union_uses_same_cell_decomposition_path():
    x, y, z = sp.symbols("x y z", real=True)
    region = sp.Or(
        sp.And(x >= 0, x <= 1, y >= 0, y <= 1, z >= 0, z <= 1),
        sp.And(x >= 2, x <= 3, y >= 0, y <= 1, z >= 0, z <= 1),
    )

    assert semialgebraic_measure(region, [x, y, z]) == 2
