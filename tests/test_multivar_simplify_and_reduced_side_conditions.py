from __future__ import annotations

import sympy as sp

from semialg.cad.lifting.lazard import lazard_evaluate, lazard_valuation
from semialg.cad.reduced import decompose_reduced_safe
from semialg.formula import Atom
from semialg.qe.complete import qe_by_complete_cad
from semialg.simplify import simp_semialg_expr
from semialg.tticad.safe import decompose_tticad_safe


def test_multivar_simpl_01():
    x, y = sp.symbols("x y", real=True)
    assert simp_semialg_expr((x > 1) & (x > 0)) == (x > 1)
    assert simp_semialg_expr(sp.Eq(x, 0) & (x**2 + y > 0)) == (sp.Eq(x, 0) & (y > 0))


def test_multivar_simpl_02():
    x, y, z = sp.symbols("x y z", real=True)
    matrix = Atom(y**2 - x, ">=") & Atom(z**2 + 1, ">")
    result = qe_by_complete_cad([z, y, x], [("exists", z)], matrix, free_variables=[x, y])
    assert result.status == "complete"
    assert result.cell_union is not None
    assert result.cell_union.cells_by_level
    assert result.formula != sp.false


def test_multivar_simpl_03():
    x, y = sp.symbols("x y", real=True)
    evaluation = lazard_evaluate((y - x) ** 2 * (y + 1), [y], [x])
    assert evaluation.valuation == (2,)
    assert sp.expand(evaluation.final_expr - (x + 1)) == 0
    assert lazard_valuation((y - x) ** 2 * (y + 1), [y], [x]) == (2,)


def test_multivar_simpl_04():
    x, y = sp.symbols("x y", real=True)
    result = decompose_reduced_safe(
        [y**2 - x, y - 1], [x, y], backend="mccallum", equational_constraints=[y**2 - x]
    )
    assert result.side_conditions is not None
    assert "well-orientedness/nullification scan" in result.side_conditions.checked_conditions
    assert result.validity is not None


def test_multivar_simpl_05():
    x, y = sp.symbols("x y", real=True)
    formula = Atom(y**2 - x, "=") | Atom(y - 1, ">")
    result = decompose_tticad_safe(formula, [x, y])
    assert result.side_conditions is not None
    assert "well-orientedness/nullification scan" in result.side_conditions.checked_conditions
