import sympy as sp

from semialg import solve_semialgebraic
from semialg.decision import _metadata


def test_formula_output_skips_structural_metadata():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    _metadata.clear_solution_metadata_cache()

    result = solve_semialgebraic(formula, (x, y), count=0, output="formula")

    assert result == formula
    cached = _metadata.collect_solution_metadata(
        formula, (x, y), request=_metadata.MetadataRequest()
    )
    assert cached["cells"] == ()
    assert cached["cylindrical_solution"] is None
    assert cached["connectivity"] is None


def test_structured_result_preserves_structural_metadata_contract():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(x >= 0, x <= 1, y >= x, y <= 1)
    _metadata.clear_solution_metadata_cache()

    result = solve_semialgebraic(formula, (x, y), count=0)

    assert result.cells
    assert result.cylindrical_solution is not None
    assert result.connectivity is None


def test_metadata_cache_reuses_cheap_analysis():
    x = sp.Symbol("x", real=True)
    formula = x**2 <= 1
    _metadata.clear_solution_metadata_cache()
    request = _metadata.MetadataRequest()

    first = _metadata.collect_solution_metadata(formula, (x,), request=request)
    cache_size = len(_metadata._METADATA_CACHE)
    second = _metadata.collect_solution_metadata(formula, (x,), request=request)

    assert first == second
    assert cache_size == 1
    assert len(_metadata._METADATA_CACHE) == cache_size
