import sympy as sp

from semialg import discretize_solution, solve_semialgebraic


def test_solution_explain_records_standard_diagnostics_for_interval_path():
    x = sp.symbols("x", real=True)
    sol = solve_semialgebraic([x**2 > 1], [x], count=0, samples="per_component")
    explanation = sol.explain()
    assert explanation["used_interval_decomposition"] is True
    assert explanation["used_domain_normalization"] is False
    assert "raw" in explanation


def test_solution_set_operations_and_predicates():
    x = sp.symbols("x", real=True)
    left = solve_semialgebraic([x >= 0, x <= 1], [x], count=0)
    right = solve_semialgebraic([x >= sp.Rational(1, 2), x <= 2], [x], count=0)
    inter = left.intersection(right)
    assert inter.contains({x: sp.Rational(3, 4)})
    assert not inter.contains({x: sp.Rational(1, 4)})
    assert inter.is_subset_of(left)
    assert left.is_disjoint_from(x < 0)
    union = left.union(right)
    assert union.contains({x: sp.Rational(3, 2)})
    diff = left.difference(right)
    assert diff.contains({x: sp.Rational(1, 4)})
    assert not diff.contains({x: sp.Rational(3, 4)})


def test_discretize_solution_for_interval_components():
    x = sp.symbols("x", real=True)
    sol = solve_semialgebraic([x**2 > 1], [x], count=0, samples="per_component")
    data = sol.discretize(bounds=[(-3, 3)])
    assert data.source == "interval-components"
    assert data.dimension == 1
    assert len(data.segments) == 2


def test_discretize_solution_for_2d_vertical_cell():
    x, y = sp.symbols("x y", real=True)
    sol = solve_semialgebraic([x >= 0, x <= 1, y >= 0, y <= 1], [x, y], count=0)
    data = discretize_solution(sol, samples_per_curve=5)
    assert data.source == "vertical-cells-2d"
    assert data.dimension == 2
    assert len(data.polygons) >= 1


def test_solve_output_diagnostics_and_plot_data_selectors():
    x = sp.symbols("x", real=True)
    explanation = solve_semialgebraic([x >= 0, x <= 1], [x], count=0, output="diagnostics")
    assert explanation["raw"]["selected_output"] == "diagnostics"
    plot_data = solve_semialgebraic([x >= 0, x <= 1], [x], count=0, output="plot_data")
    assert plot_data.source == "interval-components"
