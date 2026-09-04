import sympy as sp

from semialg import is_satisfiable
from semialg.region_analysis import (
    region_active_boundary_strata,
    region_nonsmooth_locus,
    region_singular_locus,
)


def test_singular_locus_is_restricted_to_actual_semialgebraic_boundary():
    x, y, z = sp.symbols("x y z", real=True)
    region = sp.And(sp.Eq(x * y, 0), z >= 0)

    singular = region_singular_locus(region, (x, y, z))

    assert is_satisfiable(
        sp.And(singular, sp.Eq(x, 0), sp.Eq(y, 0), sp.Eq(z, 1)),
        (x, y, z),
    )
    assert not is_satisfiable(
        sp.And(singular, sp.Eq(x, 0), sp.Eq(y, 0), z < 0),
        (x, y, z),
    )


def test_active_inequality_boundary_is_analyzed_relative_to_equality_component():
    x, y, z = sp.symbols("x y z", real=True)
    region = sp.And(sp.Eq(x, 0), y**2 - z**3 >= 0)

    singular = region_singular_locus(region, (x, y, z))

    assert sp.simplify(singular.subs({x: 0, y: 0, z: 0})) is sp.true
    assert sp.simplify(singular.subs({x: 0, y: 1, z: 1})) is sp.false


def test_disjunctive_region_reports_only_realized_active_boundary_cells():
    x, y = sp.symbols("x y", real=True)
    region = sp.Or(x >= 0, y >= 0)

    strata = region_active_boundary_strata(region, (x, y))
    x_boundary = sp.Or(
        *(s.formula for s in strata if x in s.active_residuals or -x in s.active_residuals)
    )

    assert sp.simplify(x_boundary.subs({x: 0, y: -1})) is sp.true
    assert sp.simplify(x_boundary.subs({x: 0, y: 1})) is sp.false


def test_nonsmooth_locus_uses_realized_active_sets_for_disjunction_corner():
    x, y = sp.symbols("x y", real=True)
    region = sp.Or(x >= 0, y >= 0)

    nonsmooth = region_nonsmooth_locus(region, (x, y))

    assert sp.simplify(nonsmooth.subs({x: 0, y: 0})) is sp.true
    assert sp.simplify(nonsmooth.subs({x: 0, y: -1})) is sp.false
    assert sp.simplify(nonsmooth.subs({x: -1, y: 0})) is sp.false
