import sympy as sp

from semialg import real_algebraic_feasibility, solve_real_algebraic_set


def test_ars_positive_dimensional_circle_returns_exact_witness():
    x, y = sp.symbols("x y")
    result = real_algebraic_feasibility((x**2 + y**2 - 1,), (x, y), return_result=True)
    assert result.complete is True
    assert result.satisfiable is True
    assert result.assignment is not None
    assert sp.simplify((x**2 + y**2 - 1).subs(result.assignment)) == 0
    assert result.zero_dimensional_systems


def test_ars_positive_dimensional_empty_real_set_is_certified_unsat():
    x, y = sp.symbols("x y")
    result = real_algebraic_feasibility((x**2 + y**2 + 1,), (x, y), return_result=True)
    assert result.complete is True
    assert result.satisfiable is False
    assert result.assignment is None


def test_ars_reducible_positive_dimensional_set_uses_certified_decomposition():
    x, y = sp.symbols("x y")
    point = solve_real_algebraic_set((x * y,), (x, y))
    assert point is not None
    assert sp.simplify((x * y).subs(point)) == 0


def test_ars_rejects_nonpositive_search_limits():
    x = sp.symbols("x")
    for kwargs in ({"max_generic_attempts": 0}, {"max_generic_attempts": -1}, {"max_pieces": 0}):
        try:
            real_algebraic_feasibility((x**2 - 1,), (x,), **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected ValueError for {kwargs}")


def test_ars_zero_dimensional_backend_failure_is_explicitly_incomplete():
    from semialg.algebraic.rational_univariate import RationalUnivariateError
    from semialg.real_algebraic import _solve_real_zero_dimensional

    x = sp.symbols("x")

    def failing_solver(*args, **kwargs):
        raise RationalUnivariateError("synthetic exact-backend failure")

    points, failure = _solve_real_zero_dimensional((x**2 - 1,), (x,), solver=failing_solver)
    assert points is None
    assert failure == "rur_failure"


def test_ars_uses_regular_chain_equidimensionality_for_monomial_curve():
    x, y, z = sp.symbols("x y z", real=True)
    result = real_algebraic_feasibility((x**2 - y, x * y - z), (x, y, z), return_result=True)

    assert result.complete
    assert result.satisfiable is True
    assert result.assignment is not None
    assert all(sp.simplify(poly.subs(result.assignment)) == 0 for poly in (x**2 - y, x * y - z))
    assert result.method == "ars+rur"


def test_ars_symbolic_genericity_does_not_depend_on_attempt_budget():
    x, y = sp.symbols("x y")
    result = real_algebraic_feasibility(
        (x**2 + y**2 - 1,), (x, y), max_generic_attempts=1, return_result=True
    )
    assert result.complete is True
    assert result.satisfiable is True
    assert "exact_generic_specialization" in result.notes
    assert any(note.startswith("bad_parameter_guard_degree:") for note in result.notes)
    assert all("generic_polar_search_exhausted" not in note for note in result.notes)


def test_exact_bad_parameter_avoidance_can_leave_origin():
    from semialg.real_algebraic import _rational_point_off_hypersurface

    a, b = sp.symbols("a b")
    point = _rational_point_off_hypersurface(a * (b - 1), (a, b))
    assert sp.expand((a * (b - 1)).subs({a: point[0], b: point[1]})) != 0
