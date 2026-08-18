import pytest
import sympy as sp

from semialg.formula import parse_quant_form_text
from semialg.qe.complete import CellUnion, QEDiagnostics, qe_by_complete_cad
from semialg.solve.complete import reduce_complete_text

pytestmark = pytest.mark.slow


def test_complete_qe_01():
    x = sp.Symbol("x", real=True)
    y = sp.Symbol("y", real=True)
    parsed = parse_quant_form_text(
        "exists y. y^2 - x = 0",
        symbols={"x": x, "y": y},
        variable_order=[y, x],
    )
    result = qe_by_complete_cad(parsed.vars, parsed.quantifiers, parsed.matrix)
    assert result.variables == (x, y)
    assert result.free_variables == (x,)
    assert result.quantified_variables == (y,)
    assert result.diagnostics is not None
    assert result.diagnostics.variable_reordered is True
    assert result.formula != sp.false
    assert bool(result.formula.subs(x, 1)) is True
    assert bool(result.formula.subs(x, -1)) is False


def test_complete_qe_02():
    x = sp.Symbol("x", real=True)
    y = sp.Symbol("y", real=True)
    true_result = reduce_complete_text(
        "exists x. forall y. x^2 + y^2 >= 0",
        symbols={"x": x, "y": y},
        return_result=True,
    )
    assert true_result.result == sp.true
    assert true_result.qe_result is not None
    assert true_result.qe_result.truth_value is True
    assert true_result.metadata["quantifier_blocks"] == (("exists", ("x",)), ("forall", ("y",)))

    false_result = reduce_complete_text(
        "exists x. forall y. y^2 + x = 0",
        symbols={"x": x, "y": y},
        return_result=True,
    )
    assert false_result.result == sp.false
    assert false_result.qe_result is not None
    assert false_result.qe_result.truth_value is False


def test_complete_qe_03():
    x = sp.Symbol("x", real=True)
    y = sp.Symbol("y", real=True)
    solved = reduce_complete_text(
        "exists y. y^2 - x = 0",
        symbols={"x": x, "y": y},
        variable_order=[y, x],
        return_result=True,
    )
    qe = solved.qe_result
    assert qe is not None
    assert isinstance(qe.cell_union, CellUnion)
    assert isinstance(qe.diagnostics, QEDiagnostics)
    assert qe.cell_union.variables == (x,)
    assert qe.cell_union.cell_indices == qe.satisfying_cell_indices
    assert solved.metadata["cell_union_cell_count"] == len(qe.cell_union.cells)
    assert solved.metadata["variable_reordered"] is True
