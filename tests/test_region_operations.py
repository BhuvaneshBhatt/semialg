import sympy as sp

from semialg import (
    connected_components,
    region_boundary,
    region_closure,
    region_complement,
    region_difference,
    region_dimension,
    region_interior,
    region_intersection,
    region_union,
)


def same_logic(lhs, rhs):
    return sp.simplify_logic(sp.Equivalent(lhs, rhs), form="dnf") == sp.true


def test_boolean_region_operations():
    x = sp.Symbol("x", real=True)
    assert same_logic(region_union(x < 0, x > 1), sp.Or(x < 0, x > 1))
    assert same_logic(region_intersection(x >= 0, x <= 1), sp.And(x >= 0, x <= 1))
    assert same_logic(region_difference(x >= 0, x > 1), sp.And(x >= 0, x <= 1))
    assert same_logic(region_complement(x > 0), x <= 0)


def test_interval_interior_of_closure_boundary_dimension_components():
    x = sp.Symbol("x", real=True)
    region = sp.And(x > 0, x < 1)
    assert same_logic(region_closure(region, [x]), sp.And(x >= 0, x <= 1))
    assert same_logic(region_interior(sp.And(x >= 0, x <= 1), [x]), sp.And(x > 0, x < 1))
    assert same_logic(region_boundary(region, [x]), sp.Or(sp.Eq(x, 0), sp.Eq(x, 1)))
    assert region_dimension(region, [x]) == 1
    assert len(connected_components(sp.Or(x < -1, x > 1), [x])) == 2


def test_disk_interior_of_closure_boundary_dimension():
    x, y = sp.symbols("x y", real=True)
    open_disk = x**2 + y**2 < 1
    closed_disk = x**2 + y**2 <= 1
    from semialg.decision import equivalent

    assert equivalent(region_closure(open_disk, (x, y)), closed_disk, (x, y))
    assert equivalent(region_interior(closed_disk, (x, y)), open_disk, (x, y))
    assert equivalent(region_boundary(open_disk, (x, y)), sp.Eq(x * x + y * y, 1), (x, y))
    assert region_dimension(open_disk, [x, y]) == 2
    assert region_dimension(sp.Eq(x**2 + y**2, 1), [x, y]) == 1
    assert region_dimension(sp.And(sp.Eq(x, 0), sp.Eq(y, 0)), [x, y]) == 0


def test_components_for_explicit_disjunctions_are_merged_when_touching():
    x = sp.Symbol("x", real=True)
    components = connected_components(sp.Or(sp.And(x >= 0, x <= 1), sp.And(x >= 1, x <= 2)), [x])
    assert len(components) == 1
    assert same_logic(components[0], sp.And(x >= 0, x <= 2))
