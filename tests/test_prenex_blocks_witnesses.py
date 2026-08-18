import sympy as sp

from semialg.parser import parse_quantified_formula
from semialg.qe.blocks import QuantifierBlock, norm_quant_blocks, quantifiers_to_blocks
from semialg.qe.prenex import qe_blocks, qe_text
from semialg.qe.witnesses import witness_from_qe_result
from semialg.simplify import simplify_qe_formula as canonicalize_qe_formula
from semialg.simplify.intervals import Interval1D, intervals_to_formula, merge_intervals
from semialg.validation.equivalence import sym_diff_empty


def test_prenex_blocks_01():
    x, y, z = sp.symbols("x y z", real=True)
    blocks = (
        QuantifierBlock("exists", (y,)),
        QuantifierBlock("exists", (z,)),
        QuantifierBlock("forall", (x,)),
    )
    normalized = norm_quant_blocks(blocks)
    assert normalized == (
        QuantifierBlock("exists", (y, z)),
        QuantifierBlock("forall", (x,)),
    )
    assert quantifiers_to_blocks((("exists", y), ("exists", z), ("forall", x))) == normalized


def test_prenex_blocks_02():
    y = sp.symbols("y", real=True)
    text = "exists y. y^2 < 1"
    res_text = qe_text(text)
    matrix = parse_quantified_formula(text).matrix
    res_blocks = qe_blocks(
        vars_=(y,),
        quantifier_blocks=(QuantifierBlock("exists", (y,)),),
        matrix=matrix,
    )
    assert res_text.is_sentence and res_blocks.is_sentence
    assert bool(res_text.truth_value) is True
    assert res_text.truth_value == res_blocks.truth_value


def test_prenex_blocks_03():
    res = qe_text("exists y. y^2 < 1")
    witness = witness_from_qe_result(res, {})
    assert witness is not None
    assert witness.assignment == {}


def test_prenex_blocks_04():
    x = sp.symbols("x", real=True)
    merged = merge_intervals(
        [
            Interval1D(0, 1, True, True),
            Interval1D(1, 2, False, True),
            Interval1D(3, 3, True, True),
        ]
    )
    assert len(merged) == 2
    expr = intervals_to_formula(x, merged)
    assert bool(sp.simplify(expr.subs(x, sp.Rational(3, 2))))
    assert bool(sp.simplify(expr.subs(x, 3)))
    assert not bool(sp.simplify(expr.subs(x, sp.Rational(5, 2))))


def test_prenex_blocks_05():
    x = sp.symbols("x", real=True)
    expr1 = sp.Or(sp.And(x >= 0, x <= 1), sp.And(x > 1, x <= 2))
    expr2 = sp.And(x >= 0, x <= 2)
    canon = canonicalize_qe_formula(expr1)
    assert bool(sp.simplify(canon.subs(x, 0)))
    diff = sym_diff_empty(expr1, expr2, (x,))
    assert diff.equivalent
