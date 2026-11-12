"""Large exact workloads with independently constructed mathematical answers."""

from itertools import product
from math import prod

import pytest
import sympy as sp
from algroots import NativeF4Limits

from semialg import connected_components, is_equal, solvability_conditions
from semialg.algebraic.rational_univariate import compute_rur
from semialg.algebraic_decomposition import (
    primary_decomposition,
    verify_primary_decomposition_certificate,
)

pytestmark = pytest.mark.stress

# Dense grid bases have many coefficients. Larger primes reduce accidental
# zeros in modular supports; the wider modulus permits larger rational lifts.
GRID_LIMITS = NativeF4Limits(first_prime=1000003, max_primes=256, max_modulus_bits=8192)


def _grid(variables, degree, *, multiple=False):
    """Hide a Cartesian grid behind the invertible matrix I + 11^T.

    The inverse sends grid coordinates r to x_i = r_i - sum(r)/(n+1).
    Powers 1 through degree give known local multiplicities when requested.
    """
    total = sum(variables)
    equations = tuple(
        sp.expand(
            sp.prod(
                (variable + total - root) ** (root + 1 if multiple else 1) for root in range(degree)
            )
        )
        for variable in variables
    )
    points = []
    multiplicities = []
    for roots in product(range(degree), repeat=len(variables)):
        shift = sp.Rational(sum(roots), len(variables) + 1)
        points.append(tuple(root - shift for root in roots))
        multiplicities.append(prod(root + 1 for root in roots) if multiple else 1)
    return equations, tuple(points), tuple(multiplicities)


def _check_grid_rur(dimension, degree):
    """Verify a grid RUR against its point-derived eliminant and residue identities."""
    variables = sp.symbols(f"x:{dimension}")
    source, points, _ = _grid(variables, degree)
    result = compute_rur(source, variables, domain=sp.QQ, groebner_limits=GRID_LIMITS)
    assert result.dimension == result.solution_count == degree**dimension
    defining = result.defining_polynomial
    parameter = result.parameter
    codes = tuple(
        result.separating_linear_form.subs(dict(zip(variables, point, strict=True)))
        for point in points
    )
    assert len(set(codes)) == len(points)
    expected = sp.Poly(1, parameter, domain=sp.QQ)
    for code in codes:
        expected *= sp.Poly(parameter - code, parameter, domain=sp.QQ)
    assert defining.monic() == expected
    coordinates = result.normalized_coordinate_polynomials()
    total = sum(coordinates, sp.Poly(0, parameter, domain=sp.QQ))
    for coordinate in coordinates:
        remainder = sp.Poly(1, parameter, domain=sp.QQ)
        for root in range(degree):
            remainder = (remainder * (coordinate + total - root)).rem(defining)
        assert remainder.is_zero
    separator = result.separating_linear_form.subs(
        dict(zip(variables, (p.as_expr() for p in coordinates), strict=True)), simultaneous=True
    )
    assert sp.Poly(separator - parameter, parameter, domain=sp.QQ).rem(defining).is_zero


def _check_grid_primary(dimension, degree):
    """Verify grid support, local multiplicities and the primary certificate."""
    variables = sp.symbols(f"x:{dimension}")
    source, points, multiplicities = _grid(variables, degree, multiple=True)
    result = primary_decomposition(
        source, variables, coefficient_field=sp.QQ, groebner_limits=GRID_LIMITS
    )
    assert result.complete and result.irredundant
    assert len(result.components) == degree**dimension
    assert sorted(component.degree for component in result.components) == sorted(multiplicities)
    assert (
        sum(component.degree for component in result.components)
        == (degree * (degree + 1) // 2) ** dimension
    )
    assert verify_primary_decomposition_certificate(result.certificate)
    matched = set()
    for component in result.components:
        support = [
            index
            for index, point in enumerate(points)
            if all(
                equation.subs(dict(zip(variables, point, strict=True))) == 0
                for equation in component.radical
            )
        ]
        assert len(support) == 1
        assert component.degree == multiplicities[support[0]]
        matched.add(support[0])
    assert matched == set(range(len(points)))


def _check_curved_grid(half_roots):
    """Verify the component count and one witness per sheared grid rectangle."""
    x, y = sp.symbols("x y", real=True)
    v = y - x**2
    boundary = 2 * half_roots + 1
    horizontal = sp.prod(x - root for root in range(1, boundary))
    vertical = sp.prod(v - root for root in range(1, boundary))
    formula = sp.And(x > 0, x < boundary, v > 0, v < boundary, sp.expand(horizontal * vertical) > 0)
    components = connected_components(formula, (x, y))
    assert len(components) == half_roots**2 + (half_roots + 1) ** 2
    matched = set()
    # The shear (x,v) -> (x,v+x^2) preserves connectedness. Equal-parity
    # open grid intervals are precisely the positive-product components.
    for i, j in product(range(boundary), repeat=2):
        if i % 2 != j % 2:
            continue
        px, pv = sp.Rational(2 * i + 1, 2), sp.Rational(2 * j + 1, 2)
        assignment = {x: px, y: pv + px**2}
        support = [
            index
            for index, component in enumerate(components)
            if component.subs(assignment) is sp.true
        ]
        assert len(support) == 1
        assert support[0] not in matched
        matched.add(support[0])
    assert matched == set(range(len(components)))


def _check_parameter_feasibility(degree):
    """Check a full parameter condition, including a=b=0.

    T_degree maps [-1,1] onto [-1,1]. Thus a*T_degree(x)=b has a
    solution in that interval exactly when b^2 <= a^2, even at a=0.
    """
    x, a, b = sp.symbols("x a b", real=True)
    source = sp.And(sp.Eq(a * sp.chebyshevt(degree, x) - b, 0), x >= -1, x <= 1)
    actual = solvability_conditions(source, (x,), (a, b))
    assert is_equal(actual, b**2 <= a**2, (a, b))
    assert actual.subs({a: 0, b: 0}) is sp.true
    assert actual.subs({a: 0, b: 1}) is sp.false


def test_dense_grid_rur():
    _check_grid_rur(4, 4)


def test_multiple_grid_primary():
    _check_grid_primary(3, 3)


def test_curved_grid_topology():
    _check_curved_grid(8)


def test_exceptional_parameter_fibers():
    _check_parameter_feasibility(64)
