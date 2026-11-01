"""Independent structural and parameter-topology contracts for region consumers."""

import sympy as sp

from semialg import Ball, Box, Sphere, TransformedRegion
from semialg.region_covers import bounding_cover, optimization_cover
from semialg.region_families import region_family

x, y, a = sp.symbols("x y a", real=True)


def test_cover_translation():
    image = TransformedRegion(Box(((0, 1), (2, 3))), (x + 4, y - 5), (x, y))
    assert bounding_cover(image, (x, y)) == {x: (4, 5), y: (-3, -2)}
    chart = optimization_cover(image, (x, y)).charts[0]
    point = {p: lo for p, lo, _ in chart.bounds}
    assert tuple(expr.subs(point, simultaneous=True) for expr in chart.mapping) == (4, -3)


def test_degenerate_bounds():
    assert bounding_cover(Ball((2, 3), 0), (x, y)) == {x: (2, 2), y: (3, 3)}


def test_cover_parameter_guard():
    sphere = Sphere((0, 0), a, assumptions=a >= 0)
    chart = optimization_cover(sphere, (x, y)).charts[0]
    assert sp.simplify(chart.parameter_formula().subs(a, -1)) is sp.false
    assert sp.simplify(chart.parameter_formula().subs({a: 0, chart.parameters[0]: 0})) is sp.true


def test_radial_topology_transition():
    family = region_family(Sphere((0, 0), a, assumptions=a >= 0), (x, y), (a,))
    for parameter, expected in [(-1, (-1, ())), (0, (0, (1,))), (1, (1, (1, 1)))]:
        strata = [
            s for s in family.strata if sp.simplify(s.condition.subs(a, parameter)) is sp.true
        ]
        assert len(strata) == 1
        assert (strata[0].fiber_dimension, strata[0].betti_numbers) == expected


def test_ball_topology_contract():
    family = region_family(Ball((0, 0), a, assumptions=a >= 0), (x, y), (a,))
    positive = next(s for s in family.strata if sp.simplify(s.condition.subs(a, 2)) is sp.true)
    assert positive.fiber_dimension == 2 and positive.betti_numbers == (1, 0, 0)
    assert positive.component_count == 1


def test_family_resource_limit():
    import pytest

    from semialg.errors import ResourceLimitError

    with pytest.raises(ResourceLimitError):
        region_family(sp.Eq(x * x, a), (x,), (a,), max_cells=1)


def test_empty_family_partition():
    family = region_family(sp.false, (x,), (a,))
    assert family.empty_condition is sp.true
    assert len(family.strata) == 1 and family.strata[0].fiber_dimension == -1


def test_radial_validity():
    region = Sphere((0, 0), a, assumptions=a >= 0)
    assert sp.simplify(region.as_formula((x, y)).subs({a: -1, x: 1, y: 0})) is sp.false
    family = region_family(region, (x, y), (a,))
    assert sp.simplify(family.empty_condition.subs(a, -1)) is sp.true


def test_cached_cad_budget():
    import pytest

    from semialg.cad_algorithms.decomposition import decomp_collins_complete
    from semialg.errors import ResourceLimitError

    full = decomp_collins_complete([x * x - 2], (x,))
    assert len(full.cells_by_level[1]) == 5
    with pytest.raises(ResourceLimitError):
        decomp_collins_complete([x * x - 2], (x,), max_cells=1)
    with pytest.raises(ResourceLimitError):
        decomp_collins_complete([x * x - 2], (x,), timeout=1e-20)
    assert len(decomp_collins_complete([x * x - 2], (x,)).cells_by_level[1]) == 5


def test_ellipsoid_parameter_guard():
    from semialg import Ellipsoid

    region = Ellipsoid((0, 0), [[a, 0], [0, 1]], assumptions=a > 0)
    formula = region.as_formula((x, y))
    assert sp.simplify(formula.subs({a: 0, x: 0, y: 0}, simultaneous=True)) is sp.false
    assert sp.simplify(formula.subs({a: -1, x: 0, y: 0}, simultaneous=True)) is sp.false
    assert sp.simplify(formula.subs({a: 1, x: 0, y: 0}, simultaneous=True)) is sp.true


def test_assumed_radius():
    radius = sp.Symbol("radius", positive=True)
    family = region_family(Sphere((0, 0), radius), (x, y), (radius,))
    assert len(family.strata) == 1
    assert family.strata[0].condition is sp.true
    assert family.strata[0].fiber_dimension == 1
    assert family.strata[0].sample[radius] > 0
    assert family.strata[0].betti_numbers == (1, 1)


def test_inexact_geometry_scope():
    import pytest

    source = Ball((sp.Float("0.1"), 0), sp.Float("0.2"))
    with pytest.raises(NotImplementedError):
        bounding_cover(source, (x, y))
    with pytest.raises(NotImplementedError):
        optimization_cover(source, (x, y))

    with pytest.raises(NotImplementedError):
        region_family(source, (x, y), (a,))
