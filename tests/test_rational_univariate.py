import sympy as sp

from semialg.algebraic.rational_univariate import (
    evaluate_boolean_formula_at_point,
    evaluate_relation_at_point,
    filter_rur_solutions_by_constraints,
    sign_of_algebraic_expression,
    solve_and_filter_zero_dimensional_system_with_rur,
    solve_rur_semialgebraic_system,
)


def test_algebraic_sign_detection_for_radicals():
    assert sign_of_algebraic_expression(sp.sqrt(2) - 1) == 1
    assert sign_of_algebraic_expression(1 - sp.sqrt(2)) == -1
    assert sign_of_algebraic_expression(sp.sqrt(2) ** 2 - 2) == 0


def test_evaluate_relation_at_exact_point():
    x, y = sp.symbols("x y")
    assignment = {x: sp.sqrt(2) / 2, y: sp.sqrt(2) / 2}

    assert evaluate_relation_at_point(x > 0, assignment)
    assert evaluate_relation_at_point(x <= y, assignment)
    assert evaluate_relation_at_point(sp.Eq(x**2 + y**2, 1), assignment)
    assert not evaluate_relation_at_point(y < 0, assignment)


def test_boolean_formula_filtering_keeps_positive_diagonal_point():
    x, y = sp.symbols("x y")
    solutions = solve_rur_semialgebraic_system(
        [x**2 + y**2 - 1, x - y],
        [x, y],
        x > 0,
    )

    assert solutions == ((sp.sqrt(2) / 2, sp.sqrt(2) / 2),)


def test_boolean_formula_filtering_accepts_or_and_disequality():
    x, y = sp.symbols("x y")
    formula = sp.And(sp.Or(x < 0, y < 0), sp.Ne(y, 0))
    solutions = solve_rur_semialgebraic_system(
        [x**2 + y**2 - 1, x - y],
        [x, y],
        formula,
    )

    assert solutions == ((-sp.sqrt(2) / 2, -sp.sqrt(2) / 2),)


def test_filter_candidate_points_directly():
    x, y = sp.symbols("x y")
    points = ((-1, 1), (1, 1))
    filtered = filter_rur_solutions_by_constraints(points, [x, y], [x > 0, y >= 1])

    assert filtered == ((1, 1),)


def test_structured_solution_object_exposes_assignments_and_satisfiability():
    x, y = sp.symbols("x y")
    result = solve_and_filter_zero_dimensional_system_with_rur(
        [y - x**2, y - 1],
        [x, y],
        sp.And(x < 0, y > 0),
    )

    assert result.satisfiable
    assert result.points == ((-1, 1),)
    assert result.assignments == ({x: -1, y: 1},)
    assert result.representation.defining_polynomial.degree() == 2


def test_assignment_output_for_witness_style_callers():
    x, y = sp.symbols("x y")
    solutions = solve_rur_semialgebraic_system(
        [y - x**2, y - 1],
        [x, y],
        x > 0,
        as_assignments=True,
    )

    assert solutions == ({x: 1, y: 1},)


def test_unsatisfied_constraints_return_empty_tuple():
    x, y = sp.symbols("x y")
    solutions = solve_rur_semialgebraic_system(
        [y - x**2, y - 1],
        [x, y],
        sp.And(x > 2, y > 0),
    )

    assert solutions == ()


def test_evaluate_boolean_formula_at_point_supports_not():
    x = sp.symbols("x")
    assert evaluate_boolean_formula_at_point(sp.Not(x < 0), {x: 1})
    assert not evaluate_boolean_formula_at_point(sp.Not(x > 0), {x: 1})
