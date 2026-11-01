"""Preimages contain only inputs where every map coordinate is defined and real."""

import pytest
import sympy as sp

from semialg import SemialgebraicRegion, region_preimage
from semialg.geometry_queries import semialgebraic_preimage


@pytest.mark.parametrize(
    "expression,expected",
    [
        (lambda x: 1 / x, lambda x: sp.Ne(x, 0)),
        (lambda x: sp.sqrt(x), lambda x: x >= 0),
        (lambda x: 1 / sp.sqrt(x), lambda x: x > 0),
    ],
)
def test_full_target_domain(expression, expected):
    x, y = sp.symbols("x y", real=True)
    formula = semialgebraic_preimage(expression(x), sp.true, (x,), target_variables=(y,))
    for value in (-2, -1, 0, 1, 2):
        assert bool(formula.subs(x, value)) == bool(expected(x).subs(x, value))


def test_unused_map_coordinate():
    x, u, v = sp.symbols("x u v", real=True)
    result = semialgebraic_preimage((x, 1 / x), u >= 0, (x,), target_variables=(u, v))
    assert result.subs(x, 0) is sp.false
    assert result.subs(x, 1) is sp.true
    assert result.subs(x, -1) is sp.false


def test_piecewise_pole_is_local():
    x, y = sp.symbols("x y", real=True)
    mapping = sp.Piecewise((1 / x, x > 0), (0, True))
    result = semialgebraic_preimage(mapping, sp.true, (x,), target_variables=(y,))
    assert result is sp.true


def test_tautological_target():
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_preimage(1 / x, y * y >= 0, (x,), target_variables=(y,))
    assert result == sp.Ne(x, 0)


def test_public_preimage_domain():
    x, y = sp.symbols("x y", real=True)
    whole_line = SemialgebraicRegion(sp.true, (y,))
    result = region_preimage(whole_line, 1 / x, (x,))
    assert not result.contains((0,))
    assert result.contains((sp.Rational(1, 3),))


def test_coincident_coordinate_names():
    x = sp.Symbol("x", real=True)
    result = semialgebraic_preimage(sp.sqrt(x), sp.true, (x,), target_variables=(x,))
    assert result == (x >= 0)
