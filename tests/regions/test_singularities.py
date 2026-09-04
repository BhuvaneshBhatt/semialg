from __future__ import annotations

import sympy as sp

from semialg import SemialgebraicRegion


def test_smooth_circle_has_empty_algebraic_boundary_singular_locus():
    x, y = sp.symbols("x y", real=True)
    circle = SemialgebraicRegion(sp.Eq(x**2 + y**2, 1), (x, y))
    singular = circle.singular_locus()
    assert not singular.contains((1, 0))
    assert not singular.contains((0, 1))


def test_cusp_detects_origin_and_local_dimension():
    x, y = sp.symbols("x y", real=True)
    cusp = SemialgebraicRegion(sp.Eq(y**2, x**3), (x, y))
    assert cusp.singular_locus().contains((0, 0))
    assert cusp.local_dimension((0, 0)) == 1
    assert cusp.local_dimension((1, 1)) == 1
    assert cusp.local_dimension((-1, 0)) == -1


def test_corner_is_not_misreported_as_algebraic_hypersurface_singularity():
    x, y = sp.symbols("x y", real=True)
    quadrant_box = SemialgebraicRegion(sp.And(x >= 0, x <= 1, y >= 0, y <= 1), (x, y))
    # Current API is explicitly an algebraic-boundary singularity notion.
    assert not quadrant_box.singular_locus().contains((0, 0))
    assert quadrant_box.local_dimension((0, 0)) == 2
