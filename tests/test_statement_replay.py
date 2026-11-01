"""Public representations bind to source statements, not just valid CAD geometry."""

from dataclasses import replace

import sympy as sp

from semialg import cad, replay_certificate
from semialg.cad_algorithms import decomp_collins_complete
from semialg.partial.qe import lazy_quantifier_eliminate

from ._certificate_mutation import assert_field_replacements


def _verified(result):
    return replay_certificate(result).verified is True


def test_native_replay():
    x = sp.Symbol("x", real=True)
    result = decomp_collins_complete((x,), (x,))
    assert _verified(result)
    assert not _verified(replace(result, cells_by_level={1: result.cells[:-1]}))


def test_public_statement_fields():
    x = sp.Symbol("x", real=True)
    result = cad(x > 0, (x,), strategy="collins", return_result=True)
    alternate = cad(x < 0, (x,), strategy="collins", return_result=True)
    assert_field_replacements(
        result,
        _verified,
        {
            "formula": sp.false,
            "variables": (sp.Symbol("other", real=True),),
            "cad": decomp_collins_complete((x - 1,), (x,)),
            "cell_set": alternate.cell_set,
            "status": "unknown",
            "diagnostics": {**result.diagnostics, "input_formula": x < 0},
            "function": alternate.function,
        },
    )
    root = result.function.tree[0]
    altered_function = replace(
        result.function, tree=(replace(root, truth=not root.truth), *result.function.tree[1:])
    )
    assert not _verified(replace(result, function=altered_function))


def test_free_statement_fields():
    a, x = sp.symbols("a x", real=True)
    result = lazy_quantifier_eliminate(
        (a, x), (("exists", x),), sp.Eq(x * x, a), return_result=True
    )
    tower = replace(result.tower, original_polynomials=(x * x + a,))
    truth = dict(result.truth_by_index)
    truth[(0,)] = True
    assert_field_replacements(
        result,
        _verified,
        {
            "formula": a < 0,
            "variables": (x, a),
            "free_variables": (x,),
            "quantifiers": (("forall", x),),
            "matrix": sp.Eq(x * x, -a),
            "tower": tower,
            "cells_by_level": {**result.cells_by_level, 1: result.cells_by_level[1][:-1]},
            "truth_by_index": truth,
            "status": "unknown",
        },
    )


def test_unsupported_wrapper_scope():
    x = sp.Symbol("x", real=True)
    result = cad(x > 0, (x,), operation="closure", return_result=True)
    assert replay_certificate(result).verified is None
    assert replay_certificate(cad(sp.sqrt(x) > 1, (x,), return_result=True)).verified is None
