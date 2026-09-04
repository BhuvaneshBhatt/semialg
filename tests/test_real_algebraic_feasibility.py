import sympy as sp

from semialg import real_algebraic_feasibility, solve_real_algebraic_set


def test_ars_positive_dimensional_circle_returns_exact_witness():
    x, y = sp.symbols("x y")
    result = real_algebraic_feasibility((x**2 + y**2 - 1,), (x, y))
    assert result.complete is True
    assert result.satisfiable is True
    assert result.assignment is not None
    assert sp.simplify((x**2 + y**2 - 1).subs(result.assignment)) == 0
    assert result.zero_dimensional_systems


def test_ars_positive_dimensional_empty_real_set_is_certified_unsat():
    x, y = sp.symbols("x y")
    result = real_algebraic_feasibility((x**2 + y**2 + 1,), (x, y))
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
