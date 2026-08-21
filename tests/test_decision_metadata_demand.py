import sympy as sp

from semialg import solve_semialgebraic
from semialg.decision import _metadata


def test_formula_output_skips_structural_metadata(monkeypatch):
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    calls = {"cells": 0, "cylindrical": 0}

    original_cells = _metadata._update_from_vertical_cells
    original_cyl = _metadata._update_from_cylindrical

    def count_cells(*args, **kwargs):
        calls["cells"] += 1
        return original_cells(*args, **kwargs)

    def count_cyl(*args, **kwargs):
        calls["cylindrical"] += 1
        return original_cyl(*args, **kwargs)

    monkeypatch.setattr(_metadata, "_update_from_vertical_cells", count_cells)
    monkeypatch.setattr(_metadata, "_update_from_cylindrical", count_cyl)
    _metadata.clear_solution_metadata_cache()

    result = solve_semialgebraic(formula, (x, y), count=0, output="formula")

    assert result == formula
    assert calls == {"cells": 0, "cylindrical": 0}


def test_structured_result_preserves_structural_metadata_contract():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(x >= 0, x <= 1, y >= x, y <= 1)
    _metadata.clear_solution_metadata_cache()

    result = solve_semialgebraic(formula, (x, y), count=0)

    assert result.cells
    assert result.cylindrical_solution is not None
    assert result.connectivity is None


def test_metadata_cache_reuses_cheap_analysis(monkeypatch):
    x = sp.Symbol("x", real=True)
    formula = x**2 <= 1
    calls = {"components": 0}
    original = _metadata._update_from_components

    def counted(*args, **kwargs):
        calls["components"] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(_metadata, "_update_from_components", counted)
    _metadata.clear_solution_metadata_cache()
    request = _metadata.MetadataRequest()

    first = _metadata.collect_solution_metadata(formula, (x,), request=request)
    second = _metadata.collect_solution_metadata(formula, (x,), request=request)

    assert first == second
    assert calls["components"] == 1
