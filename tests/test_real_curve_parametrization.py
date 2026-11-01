"""Real normalization coverage includes infinity and can miss isolated real points."""

from dataclasses import replace

import pytest
import sympy as sp
from algroots.curves import PlaneCurve

from semialg.curve_parametrization import verify_real_curve_parametrization

x, y, t = sp.symbols("x y t", real=True)


@pytest.mark.parametrize("f", [x * x + y * y - 1, y * y - x**3, y * y - x * x * (x + 1)])
def test_real_circle_cusp_and_two_real_node_branches(f):
    P = PlaneCurve.from_polynomial(f, (x, y)).rational_parametrization(t).parametrization
    result = verify_real_curve_parametrization(P)
    assert result.certified and result.covers_curve
    assert result.missing_condition is sp.false
    if f == x * x + y * y - 1:
        assert result.infinity_point == (0, -1)
        assert result.image_formula.subs({x: 0, y: -1}) is sp.true


def test_isolated_real_node_is_not_covered_by_real_normalization():
    C = PlaneCurve.from_polynomial(y * y - x * x * (x - 1), (x, y))
    P = C.rational_parametrization(t).parametrization
    result = verify_real_curve_parametrization(P)
    assert result.certified and not result.covers_curve
    assert result.comparison.counterexample is not None
    assert result.missing_condition.subs({x: 0, y: 0}) is sp.true
    assert result.image_formula.subs({x: 0, y: 0}) is sp.false
    restricted = verify_real_curve_parametrization(P, domain=sp.Ne(x, 0))
    assert restricted.covers_curve and restricted.certified


def test_output_renaming_and_real_extension_embedding():
    a, b = sp.symbols("a b", real=True)
    C = PlaneCurve.from_polynomial(x * x + y * y - 3, (x, y))
    P = C.rational_parametrization(t, extend_constants=True).parametrization
    result = verify_real_curve_parametrization(P, coordinates=(a, b))
    assert result.covers_curve and result.certified
    assert result.image_formula.free_symbols <= {a, b}


def test_invalid_map_complex_constants_and_capture():
    P = PlaneCurve.from_polynomial(y * y - x**3, (x, y)).rational_parametrization(t).parametrization
    with pytest.raises(ValueError, match="verified"):
        verify_real_curve_parametrization(replace(P, inverse=x))
    with pytest.raises(ValueError, match="domain"):
        verify_real_curve_parametrization(P, domain=t > 0)
    with pytest.raises(ValueError):
        verify_real_curve_parametrization(P, coordinates=(t, y))
    C = PlaneCurve.from_polynomial(x * x + y * y + 1, (x, y))
    complex_map = C.rational_parametrization(t, extend_constants=True).parametrization
    with pytest.raises(ValueError, match="real embedding"):
        verify_real_curve_parametrization(complex_map)
