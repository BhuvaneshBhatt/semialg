import sympy as sp

from semialg import IntervalComponent, solve_semialgebraic


def test_solve_semialgebraic_closed_interval_component():
    x = sp.symbols("x", real=True)
    sol = solve_semialgebraic([x**2 <= 1], [x], count=0)

    assert sol.satisfiable is True
    assert sol.dimension == 1
    assert sol.bounded is True
    assert sol.closed is True
    assert sol.compact is True
    assert sol.components == (IntervalComponent(x, -1, 1, True, True),)


def test_solve_semialgebraic_open_unbounded_components():
    x = sp.symbols("x", real=True)
    sol = solve_semialgebraic([x**2 > 1], [x], count=0)

    assert sol.satisfiable is True
    assert sol.dimension == 1
    assert sol.bounded is False
    assert sol.closed is False
    assert sol.compact is False
    assert sol.components == (
        IntervalComponent(x, -sp.oo, -1, False, False),
        IntervalComponent(x, 1, sp.oo, False, False),
    )


def test_solve_semialgebraic_punctured_interval_components():
    x = sp.symbols("x", real=True)
    sol = solve_semialgebraic([x**2 <= 1, sp.Ne(x, 0)], [x], count=0)

    assert sol.components == (
        IntervalComponent(x, -1, 0, True, False),
        IntervalComponent(x, 0, 1, False, True),
    )
    assert sol.formula == sp.Or(sp.And(x >= -1, x < 0), sp.And(x > 0, x <= 1))


def test_solve_semialgebraic_point_components_are_zero_dimensional():
    x = sp.symbols("x", real=True)
    sol = solve_semialgebraic([sp.Eq(x**2 - 1, 0)], [x], count=0)

    assert sol.dimension == 0
    assert sol.bounded is True
    assert sol.closed is True
    assert sol.compact is True
    assert sol.components == (
        IntervalComponent(x, -1, -1, True, True),
        IntervalComponent(x, 1, 1, True, True),
    )


def test_solve_semialgebraic_half_open_interval_metadata():
    x = sp.symbols("x", real=True)
    sol = solve_semialgebraic([x >= 0, x < 2], [x], count=0)

    assert sol.components == (IntervalComponent(x, 0, 2, True, False),)
    assert sol.bounded is True
    assert sol.closed is False
    assert sol.compact is False


def test_solve_semialgebraic_empty_one_dimensional_system():
    x = sp.symbols("x", real=True)
    sol = solve_semialgebraic([x**2 + 1 == 0], [x], count=0)

    assert sol.satisfiable is False
    assert sol.components == ()
    assert sol.bounded is True
    assert sol.closed is True
    assert sol.compact is True
