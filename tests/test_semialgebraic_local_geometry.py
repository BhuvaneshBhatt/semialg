import sympy as sp

import semialg

x, y = sp.symbols("x y", real=True)


def test_constraints_and_active_boundary():
    region = sp.And(sp.Eq(y, 0), x >= 0, x <= 1)
    system = semialg.polynomial_constraints(region, (x, y))
    assert len(system.clauses) == 1
    active = semialg.active_constraints(region, (0, 0), (x, y))
    assert {c.relation for c in active.constraints} == {"eq", "ge"}


def test_relative_topology_on_segment():
    region = sp.And(sp.Eq(y, 0), x >= 0, x <= 1)
    interior = semialg.relative_interior(region, (x, y))
    boundary = semialg.relative_boundary(region, (x, y))
    assert semialg.equivalent(interior, sp.And(sp.Eq(y, 0), x > 0, x < 1), (x, y))
    assert semialg.equivalent(
        boundary, sp.And(sp.Eq(y, 0), sp.Or(sp.Eq(x, 0), sp.Eq(x, 1))), (x, y)
    )


def test_semialgebraic_tangent_cone_halfline():
    cone = semialg.semialgebraic_tangent_cone(x >= 0, (0,), (x,))
    assert semialg.equivalent(cone, x >= 0, (x,))


def test_local_dimension_boundary_and_curve():
    assert semialg.local_dimension(sp.And(sp.Eq(y, x**2), x >= 0), (0, 0), (x, y)) == 1
