from __future__ import annotations

import sympy as sp

from semialg import simplify_piecewise, simplify_system
from semialg.reasoning import SimplifiedSystem
from semialg.symbolic_simplify import PiecewiseSimplificationResult


def test_simplify_piecewise_simplifies_branch_values_under_conditions():
    x = sp.symbols("x", real=True)
    expr = sp.Piecewise((sp.sqrt(x**2), x >= 0), (-x, True), evaluate=False)
    assert simplify_piecewise(expr, [x]) == sp.Piecewise((x, x >= 0), (-x, True), evaluate=False)


def test_simplify_piecewise_removes_unreachable_branch_after_coverage():
    x = sp.symbols("x", real=True)
    expr = sp.Piecewise((x, x > 0), (-x, x <= 0), (sp.Integer(999), x > 1), evaluate=False)
    result = simplify_piecewise(expr, [x], return_result=True)
    assert isinstance(result, PiecewiseSimplificationResult)
    assert result.expression == sp.Piecewise((x, x > 0), (-x, True), evaluate=False)
    assert result.removed_unreachable


def test_simplify_piecewise_assumption_can_collapse_to_value():
    x = sp.symbols("x", real=True)
    expr = sp.Piecewise((sp.sqrt(x**2), x >= 0), (-x, True), evaluate=False)
    assert simplify_piecewise(expr, [x], assumptions=x >= 0) == x


def test_simplify_system_records_simple_substitution():
    x, y = sp.symbols("x y", real=True)
    result = simplify_system([sp.Eq(y, x + 1), y >= 2], [x, y], return_result=True)
    assert isinstance(result, SimplifiedSystem)
    assert result.substitutions == {y: x + 1}
    assert result.formula != sp.false


def test_simplify_system_can_eliminate_simple_equalities():
    x, y = sp.symbols("x y", real=True)
    simplified = simplify_system([sp.Eq(y, x + 1), y >= 2], [x, y], eliminate_equalities=True)
    assert (
        simplified == (x >= 1)
        or simplified == ((x > 1) | sp.Eq(x, 1))
        or simplified == (sp.Eq(x, 1) | (x > 1))
    )


def test_simplify_system_constraints_output():
    x = sp.symbols("x", real=True)
    constraints = simplify_system([x > 1, x >= 0], [x], output="constraints")
    assert constraints == (x > 1,)
