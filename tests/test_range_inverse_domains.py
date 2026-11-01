"""Rational inverse poles and exceptional coefficient strata stay explicit."""

from itertools import product

import pytest
import sympy as sp

from semialg import function_range
from semialg._function_graph_image import (
    _try_complete_cad_image,
    _try_solved_graph_image,
)

x, y, a = sp.symbols("x y a", real=True)


def image_set(result):
    return sp.Intersection(result.formula.as_set(), sp.S.Reals)


@pytest.mark.parametrize("offset,scale,pole", [(0, 1, 0), (2, 3, -1), (-3, -2, 4)])
@pytest.mark.parametrize("method", ["qe", "cad", "auto"])
def test_reciprocal_images(offset, scale, pole, method):
    f = offset + scale / (x - pole)
    result = function_range(
        f, sp.Ne(x, pole), [x], value_symbol=y, method=method, return_result=True
    )
    expected = sp.S.Reals - sp.FiniteSet(offset)
    assert image_set(result) == expected
    assert result.formula.subs(y, offset) is sp.false
    assert result.infimum == -sp.oo and result.supremum == sp.oo
    assert result.is_interval is False and result.interval_count == 2


@pytest.mark.parametrize("f,hole,pole", [((x + 1) / (x - 1), 1, 1), ((2 * x + 3) / (x - 2), 2, 2)])
def test_fractional_linear_images(f, hole, pole):
    result = function_range(
        f, sp.Ne(x, pole), [x], value_symbol=y, method="cad", return_result=True
    )
    assert image_set(result) == sp.S.Reals - sp.FiniteSet(hole)


@pytest.mark.parametrize(
    "condition,expected",
    [
        (x > 0, sp.Interval.open(0, sp.oo)),
        (
            sp.And(x >= -1, x <= 1, sp.Ne(x, 0)),
            sp.Union(sp.Interval(-sp.oo, -1), sp.Interval(1, sp.oo)),
        ),
        (
            sp.Or(x <= -1, x >= 1),
            sp.Union(sp.Interval.Ropen(-1, 0), sp.Interval.Lopen(0, 1)),
        ),
    ],
)
def test_restricted_reciprocals(condition, expected):
    result = function_range(1 / x, condition, [x], value_symbol=y, method="cad", return_result=True)
    assert image_set(result) == expected


def test_inverse_guards():
    condition = _try_solved_graph_image(1 / x, sp.Ne(x, 0), (x,), y)
    assert condition == sp.Ne(y, 0)


def test_parameter_degree_loss():
    assert _try_solved_graph_image(a * x, sp.true, (x,), y) is None
    assert _try_solved_graph_image((a * x + 1) / (x + 1), sp.Ne(x, -1), (x,), y) is None
    nonzero = sp.Symbol("nonzero", nonzero=True, real=True)
    assert _try_solved_graph_image(nonzero * x, sp.true, (x,), y) is sp.true


def test_nonrational_inverses():
    assert _try_solved_graph_image(sp.sin(x), sp.Ne(x, 0), (x,), y) is None
    assert _try_solved_graph_image(sp.exp(x), sp.Ne(x, 0), (x,), y) is None


@pytest.mark.parametrize("method", ["qe", "cad", "auto"])
def test_finite_image(method):
    z = sp.Symbol("z", real=True)
    result = function_range(
        x,
        sp.Eq((x - 1) * (x - 3), 0),
        [x, z],
        value_symbol=y,
        method=method,
        return_result=True,
    )
    assert image_set(result) == sp.FiniteSet(1, 3)
    assert result.formula.subs(y, 2) is sp.false
    assert result.diagnostics["range_kind"] == "exact_image"


def test_bounds_are_enclosures():
    z = sp.Symbol("z", real=True)
    result = function_range(
        x,
        sp.Eq((x - 1) * (x - 3), 0),
        [x, z],
        value_symbol=y,
        method="bounds",
        return_result=True,
    )
    assert image_set(result) == sp.Interval(1, 3)
    assert result.diagnostics["range_kind"] == "interval_enclosure"
    assert result.infimum == 1 and result.supremum == 3


@pytest.mark.parametrize(
    "a,b,c,d",
    [v for v in product((-1, 0, 1), repeat=4) if v[0] * v[3] - v[1] * v[2] != 0],
)
def test_invertible_mobius(a, b, c, d):
    f = (a * x + b) / (c * x + d)
    result = function_range(
        f, sp.Ne(c * x + d, 0), [x], value_symbol=y, method="qe", return_result=True
    )
    expected = sp.S.Reals if c == 0 else sp.S.Reals - sp.FiniteSet(sp.Rational(a, c))
    assert image_set(result) == expected
    assert result.diagnostics["range_kind"] == "exact_image"


def test_complex_inverse_declines():
    assert _try_solved_graph_image(sp.I * x, sp.true, (x,), y) is None
    assert _try_solved_graph_image(x + sp.I, sp.true, (x,), y) is None


@pytest.mark.parametrize(
    "f,condition,samples",
    [
        (a * x, sp.true, [(0, sp.FiniteSet(0)), (1, sp.S.Reals), (-1, sp.S.Reals)]),
        (
            (a * x + 1) / (x + 1),
            sp.Ne(x, -1),
            [
                (0, sp.S.Reals - sp.FiniteSet(0)),
                (1, sp.FiniteSet(1)),
                (2, sp.S.Reals - sp.FiniteSet(2)),
            ],
        ),
    ],
)
def test_free_parameter_projection(f, condition, samples):
    result = _try_complete_cad_image(f, condition, (x,), y)
    assert isinstance(result, sp.Basic)
    for parameter, expected in samples:
        formula = result.subs(a, parameter)
        assert sp.Intersection(formula.as_set(), sp.S.Reals) == expected
