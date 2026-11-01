"""Exact specialized extrema retain feasibility, limits and attainment."""

from dataclasses import replace
from unittest.mock import patch

import pytest
import sympy as sp

from semialg import replay_certificate, semialgebraic_maximize, semialgebraic_minimize
from semialg._optimization_ellipsoid import linear_ellipsoid_optimum
from semialg._optimization_univariate import univariate_polynomial_optimum


@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("kind", ["min", "max"])
def test_ellipsoid_support(strict, kind):
    x, y = sp.symbols("x y", real=True)
    region = (x - 2) ** 2 + 4 * (y + 1) ** 2 < 9 if strict else (x - 2) ** 2 + 4 * (y + 1) ** 2 <= 9
    solver = semialgebraic_minimize if kind == "min" else semialgebraic_maximize
    with (
        patch(
            "semialg.optimization._certify_candidate_by_qe",
            side_effect=AssertionError("CAD decision"),
        ),
        patch(
            "semialg.optimization._certify_optimum_by_range",
            side_effect=AssertionError("image CAD"),
        ),
    ):
        result = solver(3 * x + 4 * y + 5, region, (x, y), return_result=True)
    sign = -1 if kind == "min" else 1
    assert sp.simplify(result.value - (7 + sign * 3 * sp.sqrt(13))) == 0
    assert result.attained is (not strict)
    assert result.certified
    if strict:
        assert result.points == ()
    else:
        assert sp.simplify(region.subs(result.point)) is sp.true
        assert sp.simplify(3 * result.point[x] + 4 * result.point[y] + 5 - result.value) == 0


@pytest.mark.parametrize("level", [0, 1])
def test_constant_on_ellipsoid_surface(level):
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_minimize(
        7, sp.Eq(2 * x * x + 2 * x * y + 3 * y * y, level), (x, y), return_result=True
    )
    assert result.value == 7
    assert result.attained and result.certified
    assert sp.simplify((2 * x * x + 2 * x * y + 3 * y * y - level).subs(result.point)) == 0


def test_ellipsoid_declines_extra_constraints_and_singular_forms():
    x, y, a = sp.symbols("x y a", real=True)
    assert (
        linear_ellipsoid_optimum(x, sp.And(x * x + y * y <= 1, x >= 0), (x, y), kind="min") is None
    )
    assert linear_ellipsoid_optimum(x, y * y <= 1, (x, y), kind="min") is None
    assert linear_ellipsoid_optimum(x, a * x * x + y * y <= 1, (x, y), kind="min") is None


@pytest.mark.parametrize(
    "expression,domain,kind,value,attained",
    [
        (lambda x: x * x, lambda x: sp.Ne(x, 0), "min", 0, False),
        (lambda x: x * x, lambda x: sp.Ne(x, 0), "max", sp.oo, False),
        (lambda x: -(x**3), lambda x: sp.true, "min", -sp.oo, False),
        (lambda x: x**4 - 2 * x * x, lambda x: sp.true, "min", -1, True),
        (lambda x: x, lambda x: sp.Eq(x * x, 2), "min", -sp.sqrt(2), True),
        (lambda x: sp.S(5), lambda x: x > 0, "min", 5, True),
        (lambda x: x * x, lambda x: sp.And(x > 0, x <= 2), "max", 4, True),
    ],
)
def test_polynomial_calculus(expression, domain, kind, value, attained):
    x = sp.Symbol("x", real=True)
    solver = semialgebraic_minimize if kind == "min" else semialgebraic_maximize
    with (
        patch(
            "semialg.optimization._certify_candidate_by_qe",
            side_effect=AssertionError("CAD decision"),
        ),
        patch(
            "semialg.optimization._certify_optimum_by_range",
            side_effect=AssertionError("image CAD"),
        ),
    ):
        result = solver(expression(x), domain(x), (x,), return_result=True)
    assert result.value == value
    assert result.attained is attained
    assert result.certified
    for point in result.points:
        assert sp.simplify(domain(x).subs(point)) is sp.true
        assert sp.simplify(expression(x).subs(point) - value) == 0


def test_calculus_declines_incomplete_root_coverage():
    x = sp.Symbol("x", real=True)
    with patch(
        "semialg._optimization_univariate.finite_real_roots",
        side_effect=NotImplementedError("forced"),
    ):
        assert univariate_polynomial_optimum(x * x, sp.true, x, kind="min") is None


def test_replay_rejects_infeasible_optimizer():
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_minimize(x, x * x + y * y <= 1, (x, y), return_result=True)
    assert replay_certificate(result).verified
    forged = replace(result, points=({x: -1, y: 100},))
    assert replay_certificate(forged).verified is False


def test_algebraic_ellipsoid_support():
    x, y = sp.symbols("x y", real=True)
    domain = x * x + 2 * sp.sqrt(2) * x * y + 3 * y * y <= 1
    result = semialgebraic_minimize(y, domain, (x, y), return_result=True)
    assert result.value == -1 and result.attained and result.certified
    assert result.point == {x: sp.sqrt(2), y: -1}
    assert sp.simplify(domain.subs(result.point)) is sp.true


def test_unavailable_support_sign_declines():
    x, y = sp.symbols("x y", real=True)
    with patch("semialg._optimization_ellipsoid.exact_sign", side_effect=ValueError("forced")):
        assert linear_ellipsoid_optimum(x, x * x + y * y <= 1, (x, y), kind="min") is None


def test_replay_rejects_points_for_unattained_bound():
    x = sp.Symbol("x", real=True)
    result = semialgebraic_minimize(x, x > 0, (x,), return_result=True)
    assert replay_certificate(result).verified
    assert replay_certificate(replace(result, points=({x: 0},))).verified is False


@pytest.mark.parametrize(
    "assumptions,objective,upper,kind,value,attained",
    [
        ({"positive": True}, lambda x: x, True, "min", 0, False),
        ({"negative": True}, lambda x: x, False, "max", 0, False),
        ({"nonnegative": True}, lambda x: x, True, "min", 0, True),
        ({"nonpositive": True}, lambda x: x, False, "max", 0, True),
        ({"nonzero": True, "real": True}, lambda x: x * x, True, "min", 0, False),
    ],
)
def test_optimizer_preserves_symbol_domain(assumptions, objective, upper, kind, value, attained):
    x = sp.Symbol("x", **assumptions)
    constraint = x <= 1 if upper else x >= -1
    solver = semialgebraic_minimize if kind == "min" else semialgebraic_maximize
    result = solver(objective(x), constraint, (x,), return_result=True)
    assert result.value == value and result.attained is attained and result.certified
    assert replay_certificate(result).verified
    assert set(result.variables) == {x}
    assert all(set(point) == {x} for point in result.points)


def test_replay_checks_symbol_domain():
    x = sp.Symbol("x", positive=True)
    result = semialgebraic_minimize(7, x <= 1, (x,), return_result=True)
    assert result.value == 7 and result.attained and result.certified
    assert result.point[x] > 0 and result.point[x] <= 1
    assert replay_certificate(result).verified
    assert not replay_certificate(replace(result, points=({x: -1},))).verified
