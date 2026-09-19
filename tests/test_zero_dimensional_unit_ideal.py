import sympy as sp

from semialg.solve import solve_zero_dimensional_system


def test_unit_ideal_is_unsatisfiable_not_empty_input():
    x = sp.Symbol("x")
    result = solve_zero_dimensional_system((sp.Integer(1),), variables=(x,), return_result=True)
    assert result.status == "unsat"
    assert result.points == ()
