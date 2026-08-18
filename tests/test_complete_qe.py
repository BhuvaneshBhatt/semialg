import pytest
import sympy as sp

from semialg.formula import parse_quant_form_text
from semialg.qe.complete import qe_by_complete_cad
from semialg.solve.complete import reduce_complete_text

pytestmark = pytest.mark.slow


def test_complete_qe_01():
    x = sp.Symbol("x", real=True)
    parsed = parse_quant_form_text("exists x. x^2 - 1 = 0", symbols={"x": x})
    result = qe_by_complete_cad(parsed.vars, parsed.quantifiers, parsed.matrix)
    assert result.is_sentence
    assert result.truth_value is True
    assert result.formula == sp.true


def test_complete_qe_02():
    x = sp.Symbol("x", real=True)
    parsed = parse_quant_form_text("forall x. x^2 + 1 > 0", symbols={"x": x})
    result = qe_by_complete_cad(parsed.vars, parsed.quantifiers, parsed.matrix)
    assert result.is_sentence
    assert result.truth_value is True


def test_complete_qe_03():
    x = sp.Symbol("x", real=True)
    y = sp.Symbol("y", real=True)
    result = reduce_complete_text(
        "exists y. y^2 - x = 0", symbols={"x": x, "y": y}, variable_order=[x, y]
    )
    assert result != sp.false
    assert bool(result.subs(x, 1)) is True
    assert bool(result.subs(x, -1)) is False


def test_complete_qe_04():
    x = sp.Symbol("x", real=True)
    solved = reduce_complete_text("exists x. x^2 = 1", symbols={"x": x}, return_result=True)
    assert solved.status == "complete"
    assert solved.backend == "collins-complete-qe"
    assert solved.metadata["is_sentence"] is True
