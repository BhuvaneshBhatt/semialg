"""Reduced CAD must be lifted in the same order used for QE propagation."""

import importlib
from dataclasses import replace

import pytest
import sympy as sp

from semialg.cad_algorithms.constants import (
    PROJECTION_LAZARD,
    PROJECTION_MCCALLUM,
    PROJECTION_TTICAD,
)
from semialg.formula import parse_quant_form_text

reduce_module = importlib.import_module("semialg.solve.reduce")


@pytest.mark.parametrize("backend", [PROJECTION_MCCALLUM, PROJECTION_LAZARD, PROJECTION_TTICAD])
@pytest.mark.parametrize("strategy", ["planner", "mccallum", "lazard", "tticad"])
def test_reduced_qe_normalizes_planner_order_before_lifting(monkeypatch, backend, strategy):
    a, x, y = sp.symbols("a x y", real=True)
    parsed = parse_quant_form_text(
        "forall y. exists x. x^2 - y^2 - a = 0",
        symbols={"a": a, "x": x, "y": y},
        variable_order=(x, a, y),
    )
    base = reduce_module._safe_strategy_selection(parsed, parsed.vars)
    monkeypatch.setattr(
        reduce_module,
        "_safe_strategy_selection",
        lambda *args: replace(base, backend=backend, variable_order=(x, a, y)),
    )
    # Exercise CAD even when the quadratic prepass could solve the formula.
    monkeypatch.setattr(reduce_module, "try_quadratic_vs_qe", lambda *args, **kwargs: None)
    formula, result, selection = reduce_module._reduce_reals(parsed, strategy=strategy)
    assert result.variables == (a, y, x)
    assert result.cad.tower.variables == result.variables
    assert selection.variable_order == result.variables
    assert result.quantifiers == (("forall", y), ("exists", x))
    assert result.free_variables == (a,)
    assert bool(formula.subs(a, 1))
    assert bool(formula.subs(a, 0))
    assert not bool(formula.subs(a, -1))
