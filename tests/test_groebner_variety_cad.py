from __future__ import annotations

import sympy as sp

from semialg import cad
from semialg.cad_algorithms.decomposition import decomp_groebner_variety
from semialg.cad_algorithms.projection.groebner import build_groebner_variety_projection
from semialg.formula import parse_formula
from semialg.qe.complete import qe_by_complete_cad


def test_groebner_projection_is_triangular_in_cad_order():
    x, y = sp.symbols("x y", real=True)
    formula = parse_formula(sp.And(sp.Eq(x - y, 0), sp.Eq(x**2 + y**2, 2)))

    projection = build_groebner_variety_projection(formula, (x, y))

    assert projection is not None
    assert projection.equality_context.dimension == 0
    assert projection.quotient_dimension == 2
    assert [p.as_expr() for p in projection.tower.level(1).polynomials] == [x**2 - 1]
    assert projection.tower.level(2).polynomials[0].as_expr().free_symbols == {x, y}
    assert projection.tower.metadata["variety_only"] is True


def test_groebner_lifting_follows_only_compatible_sections():
    x, y = sp.symbols("x y", real=True)
    formula = parse_formula(sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y**2 - 1, 0), sp.Eq(x * y - 1, 0)))
    projection = build_groebner_variety_projection(formula, (x, y))
    assert projection is not None

    decomposition = decomp_groebner_variety(projection)

    assert decomposition.backend == "groebner-variety"
    assert decomposition.complete is False
    assert [cell.sample_exprs for cell in decomposition.cells] == [(-1, -1), (1, 1)]
    assert decomposition.cell_count_by_level() == {1: 2, 2: 2}


def test_public_cad_uses_groebner_variety_and_reconstructs_full_points():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x - y, 0), sp.Eq(x**2 + y**2, 2), x > 0)

    result = cad(formula, (x, y), return_result=True)

    assert result.cad.backend == "groebner-variety"
    assert result.diagnostics["variety_only"] is True
    assert result.diagnostics["quotient_dimension"] == 2
    assert result.formula == sp.And(sp.Eq(x, 1), sp.Eq(y, 1))
    assert [cell.sample_exprs for cell in result.cell_set.cells] == [(1, 1)]


def test_collins_strategy_can_force_full_space_decomposition():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x - y, 0), sp.Eq(x**2 + y**2, 2), x > 0)

    result = cad(formula, (x, y), strategy="collins", return_result=True)

    assert result.cad.backend == "collins-complete"
    assert result.diagnostics["variety_only"] is False


def test_positive_dimensional_equalities_fall_back_to_collins():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x**2, 0), sp.Eq(x * y, 0), y > 0)

    result = cad(formula, (x, y), return_result=True)

    assert result.cad.backend == "collins-complete"
    assert sp.simplify_logic(sp.Equivalent(result.formula, sp.And(sp.Eq(x, 0), y > 0)))


def test_complex_only_finite_variety_returns_no_real_cells():
    (x,) = (sp.symbols("x", real=True),)
    result = cad(sp.Eq(x**2 + 1, 0, evaluate=False), (x,), return_result=True)

    assert result.cad.backend == "groebner-variety"
    assert result.cad.cells == ()
    assert result.formula is sp.false or result.formula == sp.false


def test_residual_inequality_filters_algebraic_sections_exactly():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x**2 - 2, 0), sp.Eq(y - x, 0), y < 0)

    result = cad(formula, (x, y), return_result=True)

    target = sp.And(sp.Eq(x, -sp.sqrt(2)), sp.Eq(y, -sp.sqrt(2)))
    assert result.cad.backend == "groebner-variety"
    assert sp.simplify_logic(sp.Equivalent(result.formula, target))
    assert len(result.cad.cells) == 1


def test_existential_qe_uses_finite_variety_cad_but_universal_does_not():
    x, y = sp.symbols("x y", real=True)
    existential_matrix = parse_formula(sp.And(sp.Eq(x - y, 0), sp.Eq(x**2 + y**2, 2), y > 0))
    existential = qe_by_complete_cad(
        (x, y),
        (("exists", y),),
        existential_matrix,
        free_variables=(x,),
        variable_order_strategy="preserve",
        use_presolve=False,
        return_result=True,
    )
    assert existential.cad.backend == "groebner-variety"
    assert existential.formula == sp.Eq(x, 1)

    universal_matrix = parse_formula(sp.And(sp.Eq(x - y, 0), sp.Eq(x**2 + y**2, 2)))
    universal = qe_by_complete_cad(
        (x, y),
        (("forall", y),),
        universal_matrix,
        free_variables=(x,),
        variable_order_strategy="preserve",
        use_presolve=False,
        return_result=True,
    )
    assert universal.cad.backend == "collins-complete"
    assert universal.formula is sp.false or universal.formula == sp.false


def test_inconsistent_common_equalities_short_circuit_to_empty_variety():
    x = sp.symbols("x", real=True)
    formula = sp.And(sp.Eq(x, 0), sp.Eq(x - 1, 0))

    result = cad(formula, (x,), return_result=True)

    assert result.cad.backend == "groebner-variety"
    assert result.diagnostics["equality_dimension"] == -1
    assert result.formula is sp.false or result.formula == sp.false
    assert result.cad.cells == ()
