import sympy as sp

from semialg import equivalent, is_satisfiable
from semialg.formula import parse_formula
from semialg.presolve import fourier_motzkin_eliminate
from semialg.qe.complete import qe_by_complete_cad
from semialg.qe.virtual_substitution.eliminate import try_quadratic_virtual_substitution_qe


def test_fourier_motzkin_matches_complete_qe_on_linear_projection():
    x, y = sp.symbols("x y", real=True)
    matrix = sp.And(x >= y + 1, x <= 3 - y)
    fm, removed = fourier_motzkin_eliminate(matrix, [x])
    cad = qe_by_complete_cad(
        [x, y], [("exists", x)], parse_formula(matrix), free_variables=[y], return_result=True
    ).formula
    assert removed == (x,)
    assert equivalent(fm, cad, [y])


def test_virtual_substitution_matches_complete_qe_at_discriminant_boundary():
    x, a = sp.symbols("x a", real=True)
    matrix = sp.Eq(x**2, a)
    vs = try_quadratic_virtual_substitution_qe([x, a], [("exists", x)], matrix)
    assert vs is not None
    cad = qe_by_complete_cad(
        [x, a], [("exists", x)], parse_formula(matrix), free_variables=[a], return_result=True
    ).formula
    assert equivalent(vs.formula, cad, [a])


def test_strict_and_nonstrict_touching_boundaries_are_distinguished():
    x = sp.Symbol("x", real=True)
    assert not is_satisfiable(sp.And(x > 0, x <= 0), [x])
    assert is_satisfiable(sp.And(x >= 0, x <= 0), [x])


def test_degree_drop_parameter_boundary_is_preserved():
    x, a = sp.symbols("x a", real=True)
    result = qe_by_complete_cad(
        [x, a],
        [("exists", x)],
        parse_formula(sp.Eq(a * x + 1, 0)),
        free_variables=[a],
        return_result=True,
    )
    assert equivalent(result.formula, sp.Ne(a, 0), [a])


def test_double_root_boundary_changes_strict_feasibility():
    x, a = sp.symbols("x a", real=True)
    result = qe_by_complete_cad(
        [x, a], [("exists", x)], parse_formula(x**2 < a), free_variables=[a], return_result=True
    )
    assert equivalent(result.formula, a > 0, [a])


def test_affine_coordinate_bound_fast_path_matches_complete_qe_oracle():
    from semialg.instances import coordinate_bounds

    x = sp.Symbol("x", real=True)
    for relation in (2 * x + 3 <= 0, -3 * x + 2 < 0, 5 * x - 1 >= 0, -2 * x - 7 > 0):
        bounds = coordinate_bounds(relation, (x,))
        assert bounds.complete and not bounds.inconsistent
        ((_, lower, upper),) = bounds.bounds
        ((_, lower_strict, upper_strict),) = bounds.strictness

        if upper is not sp.oo:
            violation = x >= upper if upper_strict else x > upper
            oracle = qe_by_complete_cad(
                [x],
                [("exists", x)],
                parse_formula(sp.And(relation, violation)),
                free_variables=[],
                return_result=True,
            )
            assert oracle.truth_value is False
        if lower is not -sp.oo:
            violation = x <= lower if lower_strict else x < lower
            oracle = qe_by_complete_cad(
                [x],
                [("exists", x)],
                parse_formula(sp.And(relation, violation)),
                free_variables=[],
                return_result=True,
            )
            assert oracle.truth_value is False
