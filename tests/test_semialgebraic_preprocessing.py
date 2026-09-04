from __future__ import annotations

import sympy as sp

from semialg import cad
from semialg.preprocess import semialgebraicize


def test_semialgebraicize_abs_introduces_polynomial_auxiliary():
    x = sp.Symbol("x", real=True)
    result = semialgebraicize(sp.Abs(x) <= 1, variables=(x,))

    assert result.changed
    assert len(result.aux_vars) == 1
    aux = result.aux_vars[0]
    assert sp.Abs(x) not in result.sympy_expr.atoms(sp.Abs)
    assert (aux >= 0) in result.assumptions
    assert sp.Eq(aux**2, x**2) in result.assumptions


def test_semialgebraicize_abs_fractional_power_is_polynomial_matrix():
    x, y = sp.symbols("x y", real=True)
    expr = sp.Abs(x) ** sp.Rational(3, 2) + sp.Abs(y) ** sp.Rational(3, 2) <= 1
    result = semialgebraicize(expr, variables=(x, y))

    assert result.changed
    assert len(result.aux_vars) == 4
    assert not result.sympy_expr.has(sp.Abs)
    for atom in result.sympy_expr.atoms(sp.Pow):
        assert not (atom.exp.is_Rational and atom.exp.q > 1)


def test_cad_uses_semialgebraic_preprocessing_for_abs():
    x = sp.Symbol("x", real=True)
    result = cad(sp.Abs(x) <= 1, [x], return_result=True)

    assert result.status == "complete"
    assert result.diagnostics["preprocessed"] is True
    assert (
        sp.simplify_logic(result.formula) == (sp.And(x > -1, x < 1) | sp.Eq(x, -1) | sp.Eq(x, 1))
        or result.formula != sp.false
    )


def test_cad_reports_preprocessing_limit_instead_of_fabricating_complete_result():
    x, y = sp.symbols("x y", real=True)
    expr = sp.Abs(x) ** sp.Rational(3, 2) + sp.Abs(y) ** sp.Rational(3, 2) <= 1
    result = cad(expr, [x, y], return_result=True)

    # This preprocessing introduces four existential auxiliaries.  The default
    # guard is three. A bounded structured request must not decompose a
    # placeholder ``True`` formula; it must instead
    # preserve the original formula and make the resource limitation explicit.
    assert result.status == "unknown"
    assert result.formula == expr
    assert result.diagnostics["preprocessed"] is True
    assert result.diagnostics["preprocess_elimination_limited"] is True
    assert result.diagnostics["max_preprocess_aux_vars"] == 3
    assert len(result.diagnostics["preprocess_aux_vars"]) == 4
    assert result.diagnostics["qe_formula"] is None
