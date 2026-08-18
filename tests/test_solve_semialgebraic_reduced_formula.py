import sympy as sp

from semialg import solve_semialgebraic


def test_reduced_formula_output_uses_one_dimensional_components():
    x = sp.symbols("x", real=True)

    formula = solve_semialgebraic([x**2 > 1], [x], count=0, output="reduced_formula")

    assert sp.simplify_logic(formula) == sp.Or(x < -1, x > 1)


def test_solution_conversion_methods_for_components_and_samples():
    x = sp.symbols("x", real=True)

    sol = solve_semialgebraic([x**2 <= 1, sp.Ne(x, 0)], [x], count=0, samples="per_component")

    assert len(sol.as_components()) == 2
    assert sol.as_formula() == sp.Or(sp.And(x >= -1, x < 0), sp.And(x > 0, x <= 1))
    assert sol.sample_points("per_component") == ({x: sp.Rational(-1, 2)}, {x: sp.Rational(1, 2)})
    summary = sol.to_dict()
    assert summary["component_count"] == 2
    assert summary["satisfiable"] is True


def test_piecewise_output_selector():
    x = sp.symbols("x", real=True)

    pw = solve_semialgebraic([x >= 0, x <= 1], [x], count=0, output="piecewise")

    assert isinstance(pw, sp.Piecewise)
    assert pw.args[0][0] == 1
    assert pw.args[-1] == (0, True)


def test_reduced_formula_for_supported_two_dimensional_cells():
    x, y = sp.symbols("x y", real=True)

    sol = solve_semialgebraic([x >= 0, x <= 1, y >= x, y <= 1], [x, y], count=0)
    reduced = solve_semialgebraic(
        [x >= 0, x <= 1, y >= x, y <= 1], [x, y], count=0, output="reduced_formula"
    )

    assert sol.as_cells()
    assert sp.simplify_logic(reduced) == sp.simplify_logic(sp.And(x >= 0, x <= 1, y >= x, y <= 1))
