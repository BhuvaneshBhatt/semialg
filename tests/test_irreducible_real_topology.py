from __future__ import annotations

import pytest
import sympy as sp

from semialg import connected_components, is_connected


def _assert_irreducible(poly, x, y):
    coeff, factors = sp.factor_list(poly, x, y)
    assert coeff != 0
    assert len(factors) == 1
    factor, multiplicity = factors[0]
    assert multiplicity == 1
    assert sp.expand(coeff * factor - poly) == 0


def test_smooth_irreducible_hyperbola_has_two_real_components():
    x, y = sp.symbols("x y", real=True)
    poly = x * y - 1
    _assert_irreducible(poly, x, y)

    components = connected_components(sp.Eq(poly, 0), (x, y))

    assert len(components) == 2
    assert not is_connected(sp.Eq(poly, 0), (x, y))


def test_singular_irreducible_curve_can_have_isolated_real_component():
    x, y = sp.symbols("x y", real=True)
    poly = y**2 - x**2 * (x - 1)
    _assert_irreducible(poly, x, y)

    components = connected_components(sp.Eq(poly, 0), (x, y))

    assert len(components) == 2
    assert not is_connected(sp.Eq(poly, 0), (x, y))


@pytest.mark.parametrize(
    ("parameter", "component_count"),
    [(-1, 1), (0, 1), (1, 2), (2, 2)],
)
def test_irreducible_singular_family_changes_real_component_count(parameter, component_count):
    x, y = sp.symbols("x y", real=True)
    poly = y**2 - x**2 * (x - parameter)
    _assert_irreducible(poly, x, y)

    components = connected_components(sp.Eq(poly, 0), (x, y))

    assert len(components) == component_count
    assert is_connected(sp.Eq(poly, 0), (x, y)) is (component_count == 1)
