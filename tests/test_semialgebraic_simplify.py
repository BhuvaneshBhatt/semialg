from __future__ import annotations

import sympy as sp

from semialg.formula import parse_formula
from semialg.qe.complete import qe_by_complete_cad
from semialg.simplify import Interval1D, intervals_to_formula, merge_intervals, simplify_qe_formula


def test_semialgebraic__01():
    x = sp.Symbol("x", real=True)
    intervals = (
        Interval1D(None, sp.Integer(0), False, False),
        Interval1D(sp.Integer(0), sp.Integer(0), True, True),
        Interval1D(sp.Integer(0), sp.Integer(2), False, False),
        Interval1D(sp.Integer(2), sp.Integer(2), True, True),
    )
    merged = merge_intervals(intervals)
    assert merged == (Interval1D(None, sp.Integer(2), False, True),)
    assert intervals_to_formula(x, intervals) == (x <= 2)


def test_semialgebraic__02():
    x, y = sp.symbols("x y", real=True)
    result = qe_by_complete_cad(
        [y, x], [("exists", y)], parse_formula(sp.Eq(y**2 - x, 0)), free_variables=[x]
    )
    assert result.formula == (x >= 0)
    assert result.cell_union is not None
    assert result.cell_union.formula == (x >= 0)


def test_semialgebraic__03():
    x, y = sp.symbols("x y", real=True)
    matrix = parse_formula(sp.And(sp.Eq(y**2 - x, 0), x <= 4))
    result = qe_by_complete_cad([x, y], [("exists", y)], matrix, free_variables=[x])
    assert sp.simplify_logic(result.formula ^ ((x >= 0) & (x <= 4))) == sp.false


def test_semialgebraic__04():
    x, y = sp.symbols("x y", real=True)
    expr = sp.And(sp.Eq(x, 0), x**2 + y > 0, x > -1)
    simplified = simplify_qe_formula(expr)
    assert sp.simplify_logic(simplified ^ sp.And(sp.Eq(x, 0), y > 0)) == sp.false
