from dataclasses import FrozenInstanceError

import pytest
import sympy as sp

from semialg import SemialgebraicContext, SemialgebraicRegion
from semialg.cad_algorithms.point_location import locate_cad_point
from semialg.decomposition.cylindrical import cad
from semialg.symbolic_regions import REqual


def test_semialgebraic_region_rejects_mismatched_context_formula():
    x = sp.symbols("x", real=True)
    context = SemialgebraicContext(x < 0, (x,))
    with pytest.raises(ValueError, match="context formula"):
        SemialgebraicRegion(x > 0, (x,), context=context)


def test_semialgebraic_region_rejects_mismatched_cad_formula():
    x = sp.symbols("x", real=True)
    result = cad(x < 0, (x,), output="formula", return_result=True)
    with pytest.raises(ValueError, match="CAD formula"):
        SemialgebraicRegion(x > 0, (x,), cad=result)


def test_semialgebraic_region_rejects_mismatched_cad_variables():
    x, y = sp.symbols("x y", real=True)
    result = cad(x < 0, (x,), output="formula", return_result=True)
    with pytest.raises(ValueError, match="CAD variables"):
        SemialgebraicRegion(y < 0, (y,), cad=result)


def test_context_lazy_properties_reuse_cached_value():
    x = sp.symbols("x", real=True)
    context = SemialgebraicContext(x**2 <= 1, (x,))
    first = context.parsed_formula
    second = context.parsed_formula
    assert first is second
    assert context._memo["parsed_formula"] is first


def test_context_structural_state_is_immutable_after_cache_creation():
    x = sp.symbols("x", real=True)
    context = SemialgebraicContext(x > 0, (x,))
    _ = context.polynomials
    with pytest.raises(FrozenInstanceError):
        context.formula = x < 0
    with pytest.raises(FrozenInstanceError):
        context.variables = ()


def test_empty_region_equality_evaluates_to_true():
    relation = REqual()
    assert relation.as_formula() is sp.true
    assert relation.evaluate() is True


def test_point_location_missing_mapping_coordinate_is_value_error():
    x, y = sp.symbols("x y", real=True)
    region = SemialgebraicRegion(sp.And(x >= 0, y >= 0), (x, y)).as_cad_region()
    with pytest.raises(ValueError, match="missing coordinate"):
        locate_cad_point(region, {x: 1})


def test_cad_function_tree_index_is_reused_and_read_only():
    x = sp.symbols("x", real=True)
    function = cad(x**2 <= 1, (x,), output="formula", return_result=True).as_function()
    first = function.tree_by_index()
    second = function.tree_by_index()
    assert first is second
    with pytest.raises(TypeError):
        first[(999,)] = next(iter(first.values()))


def test_quantified_region_validates_supplied_cad_before_reuse():
    from semialg.quantifiers import Exists

    x, z = sp.symbols("x z", real=True)
    quantified = Exists((z,), sp.Eq(z**2, x))
    wrong = cad(x < 0, (x,), output="formula", return_result=True)
    region = SemialgebraicRegion(quantified, (x,), cad=wrong)
    with pytest.raises(ValueError, match="CAD formula"):
        _ = region.cad
