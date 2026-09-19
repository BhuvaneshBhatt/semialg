from __future__ import annotations

import ast
from pathlib import Path

import sympy as sp

from semialg.decision._rur import finite_points_formula, try_rur_formula

ROOT = Path(__file__).resolve().parents[1]


def test_predicate_layer_does_not_import_decision_api() -> None:
    path = ROOT / "src" / "semialg" / "decision" / "_predicates.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imports = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    }
    assert "api" not in imports


def test_finite_rur_support_preserves_exact_formula_contract() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = finite_points_formula(({x: 1, y: 2}, {x: 3, y: 4}), (x, y))
    assert formula == sp.Or(sp.And(sp.Eq(x, 1), sp.Eq(y, 2)), sp.And(sp.Eq(x, 3), sp.Eq(y, 4)))

    result = try_rur_formula(sp.And(sp.Eq(x**2, 1), sp.Eq(y, x)), (x, y))
    assert result is not None
    assert result.partial is False
    assert {tuple(point[var] for var in (x, y)) for point in result.assignments} == {
        (-1, -1),
        (1, 1),
    }
