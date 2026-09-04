import sympy as sp

from semialg import region_boundary, region_closure, region_interior


def _truth(expr, symbol, value):
    return bool(expr.subs({symbol: value}))


def test_touching_closed_intervals_have_no_internal_boundary_seam():
    x = sp.Symbol("x", real=True)
    region = sp.Or(sp.And(x >= 0, x <= 1), sp.And(x >= 1, x <= 2))

    interior = region_interior(region, [x])
    boundary = region_boundary(region, [x])
    closure = region_closure(region, [x])

    assert _truth(interior, x, 1)
    assert _truth(interior, x, sp.Rational(1, 2))
    assert not _truth(interior, x, 0)
    assert not _truth(boundary, x, 1)
    assert _truth(boundary, x, 0)
    assert _truth(boundary, x, 2)
    assert _truth(closure, x, 1)


def test_boolean_region_hole_boundary_is_cad_semantic():
    x = sp.Symbol("x", real=True)
    region = sp.And(x >= -2, x <= 2, sp.Not(sp.And(x > -1, x < 1)))

    interior = region_interior(region, [x])
    boundary = region_boundary(region, [x])

    assert _truth(interior, x, sp.Rational(3, 2))
    assert not _truth(interior, x, 1)
    assert _truth(boundary, x, -2)
    assert _truth(boundary, x, -1)
    assert _truth(boundary, x, 1)
    assert _truth(boundary, x, 2)


def test_touching_rectangles_have_no_internal_vertical_boundary_seam():
    x, y = sp.symbols("x y", real=True)
    region = sp.Or(
        sp.And(x >= 0, x <= 1, y >= 0, y <= 1),
        sp.And(x >= 1, x <= 2, y >= 0, y <= 1),
    )

    interior = region_interior(region, [x, y])
    boundary = region_boundary(region, [x, y])

    seam_point = {x: 1, y: sp.Rational(1, 2)}
    outer_point = {x: 0, y: sp.Rational(1, 2)}
    assert bool(interior.subs(seam_point))
    assert not bool(boundary.subs(seam_point))
    assert bool(boundary.subs(outer_point))
