"""Optional downstream contract for the Semialg surface consumed by asymptotic.

This test deliberately does not make asymptotic a Semialg dependency. Release
validation can install the downstream package and select the ``downstream``
marker to exercise the real integration.
"""

import pytest
import sympy as sp

pytestmark = [pytest.mark.downstream, pytest.mark.integration]


def test_semialg_limit_support_surface_is_usable_with_downstream_asymptotic():
    pytest.importorskip("asymptotic")
    from semialg import (
        angular_map_image,
        correlated_map_image,
        local_algebraic_strata,
        point_in_closure,
    )

    x, y, u, v = sp.symbols("x y u v", real=True)
    assert point_in_closure(sp.And(x**2 + y**2 < 1, sp.Ne(x**2 + y**2, 0)), (0, 0), (x, y))
    image = angular_map_image((x, y), (x, y), image_variables=(u, v))
    assert sp.simplify(image.formula.subs({u: 1, v: 1})) is sp.false
    curve = correlated_map_image((x, x**2), sp.And(x >= -1, x <= 1), (x,), image_variables=(u, v))
    assert sp.simplify(curve.formula.subs({u: sp.Rational(1, 2), v: sp.Rational(1, 4)})) is sp.true
    strata = local_algebraic_strata((x * y,), {x: 0, y: 0}, (x, y))
    assert strata.complete and len(strata.branch_geometry.branches) == 2
