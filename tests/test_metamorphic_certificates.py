from __future__ import annotations

import sympy as sp

from semialg import convexity_certificate, is_equal
from semialg.convexity import polynomial_convexity_certificate
from semialg.instances import coordinate_bounds
from semialg.simplify.equality import simplify_equalities


def test_affine_equality_simplification_is_invariant_under_nonzero_constant_scaling():
    x, y = sp.symbols("x y", real=True)
    base = sp.And(sp.Eq(2 * x + y, 0), y >= -2, y <= 2)
    scaled = sp.And(sp.Eq(-6 * x - 3 * y, 0), y >= -2, y <= 2)
    assert is_equal(simplify_equalities(base), simplify_equalities(scaled), [x, y])


def test_coordinate_bounds_are_invariant_under_conjunct_reordering_and_redundancy():
    x = sp.Symbol("x", real=True)
    first = coordinate_bounds(sp.And(x > -2, x <= 3, x <= 5), (x,))
    second = coordinate_bounds(sp.And(x <= 5, x <= 3, x > -2), (x,))
    assert first == second
    assert first.bounds == ((x, -2, 3),)
    assert first.strictness == ((x, True, False),)


def test_constant_hessian_inertia_certificate_replays_on_diagonal_quadratic():
    x, y, z = sp.symbols("x y z", real=True)
    expr = 2 * x**2 + 3 * y**2
    cert = polynomial_convexity_certificate(expr, (x, y, z))
    hessian = sp.hessian(expr, (x, y, z))

    assert cert.outcome is True
    assert cert.hessian_inertia == (2, 0, 1)
    assert hessian == sp.diag(4, 6, 0)


def test_segment_counterexample_certificate_replays_from_returned_witness():
    x, y = sp.symbols("x y", real=True)
    horizontal = sp.And(x >= 0, x <= 2, y >= 0, y <= 1)
    vertical = sp.And(x >= 0, x <= 1, y >= 0, y <= 2)
    region = sp.Or(horizontal, vertical)
    cert = convexity_certificate(region, (x, y))

    assert cert.outcome is False
    assert cert.method == "exact-segment-counterexample"
    left = cert.details["left"]
    right = cert.details["right"]
    midpoint = cert.witness
    assert bool(region.subs(left))
    assert bool(region.subs(right))
    assert not bool(region.subs(midpoint))
    for variable in (x, y):
        assert sp.simplify(midpoint[variable] - (left[variable] + right[variable]) / 2) == 0


def test_rur_points_replay_against_original_polynomial_system():
    from semialg.algebraic.rational_univariate import solve_zero_dimensional_system_with_rur

    x, y = sp.symbols("x y", real=True)
    system = (x**2 + y**2 - 1, x - y)
    points = solve_zero_dimensional_system_with_rur(system, (x, y), real=True)
    assert points
    for point in points:
        assignment = dict(zip((x, y), point, strict=True))
        assert all(sp.simplify(poly.subs(assignment)) == 0 for poly in system)


def test_root_isolation_certificate_replays_by_exact_root_count():
    from semialg.algebraic.roots import isolate_real_roots

    x = sp.Symbol("x", real=True)
    poly = sp.Poly((x**2 - 2) * (x - 3), x)
    roots = isolate_real_roots(poly)
    assert len(roots) == 3
    for root in roots:
        value = root.as_expr()
        assert sp.simplify(poly.as_expr().subs(x, value)) == 0
        left, right = root.interval.left, root.interval.right
        assert int(poly.count_roots(left, right)) == 1


def test_optimization_certificate_replays_feasibility_value_and_no_better_point():
    from semialg import is_satisfiable, semialgebraic_minimize

    x, y = sp.symbols("x y", real=True)
    objective = x**2 + y**2
    constraints = sp.And(x + y >= 1, x >= 0, y >= 0)
    result = semialgebraic_minimize(objective, constraints, (x, y), return_result=True)
    assert result.certified is True
    assert result.point is not None
    assert bool(constraints.subs(result.point))
    assert sp.simplify(objective.subs(result.point) - result.value) == 0
    assert not is_satisfiable(sp.And(constraints, objective < result.value), (x, y))
