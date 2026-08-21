from __future__ import annotations

import pytest
import sympy as sp

import semialg.symbolic_simplify as symbolic_simplify


def test_symbolic_simplifier_does_not_swallow_programming_errors(monkeypatch):
    x = sp.Symbol("x", real=True)

    def broken_backend(*args, **kwargs):
        raise AssertionError("programming defect")

    monkeypatch.setattr(sp, "simplify_logic", broken_backend)
    with pytest.raises(AssertionError, match="programming defect"):
        symbolic_simplify._safe_simplify_logic(sp.Or(x > 0, x <= 0))
