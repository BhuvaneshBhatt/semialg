import sympy as sp

from semialg.formula import parse_formula, parse_quant_form_text
from semialg.partial.qe import lazy_find_inst_form, lazy_resolve_formula
from semialg.solve.find_instance import find_instance_text
from semialg.solve.reduce import reduce_text
from semialg.solve.resolve import resolve_text


def test_lazy_partial_01():
    parsed = parse_quant_form_text("exists x. x**2 - 1 = 0")
    result = lazy_resolve_formula(parsed.vars, parsed.quantifiers, parsed.matrix)
    assert result.truth_value is True
    assert result.stats.stopped_early is True
    assert result.stats.evaluated_leaf_cells < result.stats.visited_cells_by_level[1]
    assert result.witness is not None


def test_lazy_partial_02():
    parsed = parse_quant_form_text("forall x. x**2 - 1 = 0")
    result = lazy_resolve_formula(parsed.vars, parsed.quantifiers, parsed.matrix)
    assert result.truth_value is False
    assert result.stats.stopped_early is True
    assert result.counterexample is not None


def test_lazy_partial_03():
    x = sp.Symbol("x")
    matrix = parse_formula(sp.Eq(x**2 - 1, 0))
    result = lazy_find_inst_form((x,), matrix)
    assert result.found is True
    assert result.instance is not None
    assert sp.simplify(result.instance[x] ** 2 - 1) == 0


def test_lazy_partial_04():
    result = resolve_text("exists x. x**2 - 1 = 0", return_result=True)
    assert result.result is True
    assert result.method == "partial_cad_resolve"
    assert result.metadata["stats"].stopped_early is True


def test_lazy_partial_05():
    result = find_instance_text("x**2 - 1 = 0", return_result=True)
    assert result.method == "partial_cad_instance"
    assert result.result is not None


def test_lazy_partial_06():
    result = reduce_text("exists y. y**2 - x = 0", strategy="mccallum", return_result=True)
    assert result.result == (sp.Symbol("x") >= 0)
    selection = result.metadata["strategy_selection"]
    assert "mccallum" in str(selection.backend)
