import pytest
import sympy as sp

from semialg import Ball, Ellipsoid, Triangle
from semialg.geometry_validity import normalize_geometry_assumptions, validity_status


def test_ball_exposes_symbolic_radius_condition():
    r = sp.Symbol("r", real=True)
    ball = Ball((0, 0), r)
    assert ball.conditions == (r >= 0,)
    assert ball.is_valid() is None
    assert ball.is_valid(r >= 0) is True
    assert ball.is_valid(r < 0) is False


def test_ball_accepts_sympy_positive_predicate_assumption():
    r = sp.Symbol("r", real=True)
    ball = Ball((0, 0), r, assumptions=sp.Q.positive(r))
    assert ball.is_valid(sp.Q.positive(r)) is True
    assert normalize_geometry_assumptions(sp.Q.positive(r)) == (r > 0,)


def test_constructor_rejects_assumption_contradicting_symbolic_radius():
    r = sp.Symbol("r", real=True)
    with pytest.raises(ValueError, match="validity conditions contradict"):
        Ball((0, 0), r, assumptions=r < 0)


def test_existing_construction_conditions_feed_common_protocol():
    a, b, c = sp.symbols("a b c", positive=True)
    triangle = Triangle.from_sides(a, b, c)
    assert triangle.conditions == triangle.construction_conditions
    assert triangle.is_valid(a + b > c) is None
    assumptions = sp.And(a + b > c, a + c > b, b + c > a)
    assert triangle.is_valid(assumptions) is True


def test_ellipsoid_conditions_feed_common_protocol():
    a, b = sp.symbols("a b", real=True)
    ellipsoid = Ellipsoid((0, 0), sp.diag(a, b))
    assert ellipsoid.conditions == ellipsoid.construction_conditions
    assert ellipsoid.is_valid(sp.And(a > 0, b > 0)) is True
    assert ellipsoid.is_valid(a < 0) is False


def test_validity_status_of_unconditional_geometry():
    assert validity_status(()) is True


def test_ellipsoid_constructor_uses_shared_assumptions_contract():
    a, b = sp.symbols("a b", real=True)
    ellipsoid = Ellipsoid((0, 0), sp.diag(a, b), assumptions=sp.And(a > 0, b > 0))
    assert ellipsoid.is_valid(sp.And(a > 0, b > 0)) is True
    with pytest.raises(ValueError, match="validity conditions contradict"):
        Ellipsoid((0, 0), sp.diag(a, b), assumptions=a < 0)


def test_ellipsoid_preserves_unresolved_matrix_conditions():
    a, b, c = sp.symbols("a b c", real=True)
    ellipsoid = Ellipsoid((0, 0), ((a, b), (c, a)))
    assert (
        sp.Eq(b, c) in ellipsoid.construction_conditions
        or sp.Eq(c, b) in ellipsoid.construction_conditions
    )
    assert len(ellipsoid.construction_conditions) == 3


def test_ellipsoid_rejects_certified_nonsymmetric_numeric_matrix():
    with pytest.raises(ValueError, match="symmetric"):
        Ellipsoid((0, 0), ((2, 1), (0, 2)))
