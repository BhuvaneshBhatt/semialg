import sympy as sp

from semialg.formula import parse_formula
from semialg.heuristics import suggest_cad_variable_order
from semialg.presolve import fourier_motzkin_eliminate, presolve_semialgebraic
from semialg.qe.complete import qe_by_complete_cad


def test_presolve_affine_existential_substitution_reduces_qe_dimension():
    x, y = sp.symbols("x y", real=True)
    result = qe_by_complete_cad(
        [x, y],
        [("exists", y)],
        parse_formula(sp.And(sp.Eq(y, x + 1), y > 0)),
        free_variables=[x],
        return_result=True,
    )
    assert result.formula == (x > -1)
    assert y not in result.variables
    assert any("presolve eliminated" in note for note in result.diagnostics.notes)


def test_presolve_does_not_cross_universal_dependency():
    x, y = sp.symbols("x y", real=True)
    result = qe_by_complete_cad(
        [x, y],
        [("exists", x), ("forall", y)],
        parse_formula(sp.Eq(x, y)),
        return_result=True,
    )
    assert result.truth_value is False
    assert result.quantifiers == (("exists", x), ("forall", y))


def test_parameter_dependent_affine_coefficient_is_not_divided_away():
    x, a = sp.symbols("x a", real=True)
    result = presolve_semialgebraic(sp.Eq(a * x, 1), [a, x], eliminate=[x])
    assert result.substitutions == ()
    assert x in result.variables


def test_fourier_motzkin_exact_linear_elimination():
    x, y = sp.symbols("x y", real=True)
    reduced = fourier_motzkin_eliminate(sp.And(y > x, y < 2), [y])
    assert reduced is not None
    formula, removed = reduced
    assert removed == (y,)
    assert sp.simplify_logic(sp.Equivalent(formula, x < 2)) is sp.true


def test_presolve_reports_independent_variable_blocks():
    x, y, z = sp.symbols("x y z", real=True)
    result = presolve_semialgebraic(sp.And(x**2 <= 1, y**2 + z**2 <= 1), [x, y, z])
    assert result.variable_blocks == ((x,), (y, z))


def test_projection_scored_variable_order_is_a_permutation():
    x, y, z = sp.symbols("x y z", real=True)
    score = suggest_cad_variable_order(
        [z**4 + x * z + y, y**2 - x], [x, y, z], strategy="projection"
    )
    assert set(score.order) == {x, y, z}
    assert score.strategy == "projection"
    assert score.projection_polynomials is not None
    assert score.sum_total_degree is not None
