"""Topology rules are checked against exact sets and the forced CAD path."""

import itertools
from unittest.mock import patch

import pytest
import sympy as sp

from semialg import region_boundary, region_closure, region_interior
from semialg.decision import equivalent
from semialg.reasoning_regions import region_bounded
from semialg.regions.structural_topology import structural_topology


@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("operator", [sp.Ge, sp.Gt, sp.Le, sp.Lt])
def test_semidefinite_zero_level(sign, operator):
    x, y = sp.symbols("x y", real=True)
    polynomial = sign * sp.expand((x - y) ** 2)
    relation = operator(polynomial, 0, evaluate=False)
    for operation in (region_closure, region_interior, region_boundary):
        result = operation(relation, (x, y))
        oracle = operation(relation, (x, y), strategy="cad")
        assert equivalent(result, oracle, (x, y))


@pytest.mark.parametrize(
    "polynomial",
    [
        lambda x, y: x * x + y,
        lambda x, y: x * x - y * y,
        lambda x, y: 2 * x * x + 2 * x * y + 3 * y * y - 1,
    ],
)
def test_quadratic_sign_changing_boundary(polynomial):
    x, y = sp.symbols("x y", real=True)
    expression = polynomial(x, y)
    result = region_boundary(expression <= 0, (x, y))
    assert equivalent(result, sp.Eq(expression, 0), (x, y))
    for assignment in ({x: 0, y: 0}, {x: 1, y: 1}, {x: -1, y: 0}):
        assert bool(result.subs(assignment)) is bool(sp.Eq(expression, 0).subs(assignment))


def test_nonnegative_quartic_has_no_boundary():
    x, y = sp.symbols("x y", real=True)
    nonnegative = sp.expand((x - y) ** 4)
    assert region_boundary(sp.Ge(nonnegative, 0, evaluate=False), (x, y)) is sp.false
    assert region_interior(sp.Ge(nonnegative, 0, evaluate=False), (x, y)) is sp.true


def test_cartesian_holes_and_thin_factors():
    x, y = sp.symbols("x y", real=True)
    region = sp.And(x > 0, x < 2, sp.Ne(x, 1), sp.Eq(y, 0))
    with patch(
        "semialg.regions.operations._cad_topology_formula",
        side_effect=AssertionError("CAD topology"),
    ):
        closure = region_closure(region, (x, y))
        interior = region_interior(region, (x, y))
        boundary = region_boundary(region, (x, y))
    assert interior is sp.false
    expected = sp.And(x >= 0, x <= 2, sp.Eq(y, 0))
    assert equivalent(closure, expected, (x, y))
    assert equivalent(boundary, expected, (x, y))


def test_cartesian_boundary_faces():
    x, y, z = sp.symbols("x y z", real=True)
    region = sp.And(x > 0, x <= 1, y >= -1, y < 2, z >= 0, z <= 3)
    with patch(
        "semialg.regions.operations._cad_topology_formula",
        side_effect=AssertionError("CAD topology"),
    ):
        boundary = region_boundary(region, (x, y, z))
    for values in itertools.product((-1, 0, sp.Rational(1, 2), 1, 2, 3, 4), repeat=3):
        a, b, c = values
        expected = (
            0 <= a <= 1
            and -1 <= b <= 2
            and 0 <= c <= 3
            and (a in (0, 1) or b in (-1, 2) or c in (0, 3))
        )
        assert bool(boundary.subs(dict(zip((x, y, z), values, strict=True)))) is bool(expected)


def test_parameter_geometry_declines():
    x, y, a = sp.symbols("x y a", real=True)
    assert structural_topology(a * x * x + y * y <= 0, (x, y), "boundary") is None


def test_coercive_union_boundedness():
    x, y = sp.symbols("x y", real=True)
    union = sp.Or((x - 3) ** 2 + 2 * y * y <= 1, (x + 3) ** 2 + 2 * y * y <= 1)
    with patch(
        "semialg.reasoning_regions._unbounded_by_qe", side_effect=AssertionError("boundedness QE")
    ):
        assert region_bounded(union, (x, y)) is True
    assert region_bounded(sp.And(y >= 0, y <= 1), (x, y)) is False
