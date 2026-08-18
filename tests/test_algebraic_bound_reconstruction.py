from __future__ import annotations

import sympy as sp

from semialg.reconstruct.radicals import fiber_root_candidates, fiber_root_expr
from semialg.reconstruct.root_functions import root_of


def test_quadratic_fiber_roots_use_radicals_for_disk_boundaries():
    x, y = sp.symbols("x y", real=True)
    low, high = fiber_root_candidates(x**2 + y**2 - 1, y)
    assert sp.simplify(low + sp.sqrt(1 - x**2)) == 0
    assert sp.simplify(high - sp.sqrt(1 - x**2)) == 0


def test_quadratic_fiber_roots_support_sphere_z_bounds():
    x, y, z = sp.symbols("x y z", real=True)
    low, high = fiber_root_candidates(x**2 + y**2 + z**2 - 1, z)
    assert sp.simplify(low + sp.sqrt(1 - x**2 - y**2)) == 0
    assert sp.simplify(high - sp.sqrt(1 - x**2 - y**2)) == 0


def test_linear_fiber_root_reconstruction():
    x, y = sp.symbols("x y", real=True)
    (root,) = fiber_root_candidates(y - x**2, y)
    assert sp.simplify(root - x**2) == 0


def test_higher_degree_falls_back_to_root_of():
    x, y = sp.symbols("x y", real=True)
    expr = fiber_root_expr(y**5 + x * y + 1, y, 2)
    assert expr == root_of(x * y + y**5 + 1, y, 2)
