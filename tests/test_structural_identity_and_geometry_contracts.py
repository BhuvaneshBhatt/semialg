import pytest
import sympy as sp

import semialg
from semialg.cache_utils import BoundedLRU
from semialg.parameters import SolvabilityConditionsResult
from semialg.sampling import _dedupe_points
from semialg.symbolic_regions import (
    RegionElement,
    RSubset,
    SemialgebraicRegion,
    region_element_conditions,
    region_relation_conditions,
)
from semialg.symbolic_simplify import simplify_boole


def _same_name_symbols():
    return sp.Symbol("x"), sp.Symbol("x", real=True), sp.Symbol("x", integer=True)


@pytest.mark.parametrize("left_index,right_index", [(0, 1), (0, 2), (1, 2)])
def test_boolean_dedup_preserves_same_name_symbol_identity(left_index, right_index):
    symbols = _same_name_symbols()
    left, right = symbols[left_index], symbols[right_index]
    expr = sp.And(left > 0, right > 0, evaluate=False)
    simplified = simplify_boole(expr, (left, right), semantic=False)
    assert left in simplified.free_symbols
    assert right in simplified.free_symbols


@pytest.mark.parametrize("left_index,right_index", [(0, 1), (0, 2), (1, 2)])
def test_point_dedup_preserves_same_name_symbol_identity(left_index, right_index):
    symbols = _same_name_symbols()
    left, right = symbols[left_index], symbols[right_index]
    points = _dedupe_points(({left: 1}, {right: 1}))
    assert len(points) == 2
    assert tuple(points[0]) == (left,)
    assert tuple(points[1]) == (right,)


def test_region_element_conditions_recurse_through_implies_and_equivalent():
    x = sp.Symbol("x", real=True)
    region = SemialgebraicRegion(x >= 0, (x,))
    atom = RegionElement((x,), region)

    implication = sp.Implies(atom, x <= 2, evaluate=False)
    lowered_implication = region_element_conditions(implication)
    assert not lowered_implication.has(RegionElement)
    assert semialg.equivalent(lowered_implication, sp.Implies(x >= 0, x <= 2), (x,))

    equivalent = sp.Equivalent(atom, x >= 0, evaluate=False)
    lowered_equivalent = region_element_conditions(equivalent)
    assert not lowered_equivalent.has(RegionElement)
    assert semialg.is_tautology(lowered_equivalent, (x,))


def test_region_relation_conditions_recurse_through_implies_and_equivalent():
    x = sp.Symbol("x", real=True)
    left = SemialgebraicRegion(x >= 0, (x,))
    right = SemialgebraicRegion(x >= -1, (x,))
    subset = RSubset(left, right)

    lowered = region_relation_conditions(sp.Implies(subset, x <= 5, evaluate=False))
    assert not lowered.has(RSubset)
    lowered_eqv = region_relation_conditions(sp.Equivalent(subset, subset, evaluate=False))
    assert not lowered_eqv.has(RSubset)


def test_direct_region_membership_formula_is_structural_by_default():
    x, a = sp.symbols("x a", real=True)
    region = SemialgebraicRegion(sp.And(x >= 0, x <= a), (x,))
    atom = RegionElement((x,), region)
    assert atom.as_formula() == atom.as_formula(eliminate=False)


def test_bounded_lru_treats_cached_none_as_a_hit_for_recency():
    cache = BoundedLRU[object](2, "test-none-sentinel")
    cache.put("a", None)
    cache.put("b", 1)
    assert cache.get("a") is None
    cache.put("c", 2)
    assert tuple(cache._data) == ("a", "c")


def test_conditional_solvability_has_no_ambiguous_truth_value():
    a = sp.Symbol("a", real=True)
    result = SolvabilityConditionsResult(a >= 0, a >= 0, (), (a,))
    assert result.is_conditional
    assert not result.is_never_solvable
    assert not result.is_unconditionally_solvable
    with pytest.raises(TypeError, match="parameter-dependent solvability"):
        bool(result)


def test_unconditional_and_impossible_solvability_have_boolean_truth():
    assert bool(SolvabilityConditionsResult(sp.true, sp.true, (), ())) is True
    assert bool(SolvabilityConditionsResult(sp.false, sp.false, (), ())) is False


def test_singular_locus_uses_joint_jacobian_rank_not_individual_gradients():
    x, y, z = sp.symbols("x y z", real=True)
    region = sp.And(sp.Eq(z, 0), sp.Eq(z - x * y, 0))
    singular = semialg.region_singular_locus(region, (x, y, z))
    assert semialg.is_satisfiable(
        sp.And(singular, sp.Eq(x, 0), sp.Eq(y, 0), sp.Eq(z, 0)), (x, y, z)
    )
    assert not semialg.is_satisfiable(
        sp.And(singular, sp.Eq(x, 1), sp.Eq(y, 0), sp.Eq(z, 0)), (x, y, z)
    )


def test_active_boundary_strata_recover_only_realized_active_sets():
    x, y = sp.symbols("x y", real=True)
    region = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    strata = semialg.region_active_boundary_strata(region, (x, y))
    assert sum(s.active_count == 1 for s in strata) == 4
    assert sum(s.active_count == 2 for s in strata) == 4
    assert len(strata) == 8


def test_active_boundary_strata_scale_with_realized_box_strata_not_all_subsets():
    x, y, z = sp.symbols("x y z", real=True)
    region = sp.And(x >= 0, x <= 1, y >= 0, y <= 1, z >= 0, z <= 1)
    strata = semialg.region_active_boundary_strata(region, (x, y, z))
    counts = {
        active_count: sum(s.active_count == active_count for s in strata)
        for active_count in (1, 2, 3)
    }
    assert counts == {1: 6, 2: 12, 3: 8}
    assert len(strata) == 26


def test_symbolic_region_api_uses_pythonic_condition_names_only():
    import semialg.symbolic_regions as regions

    assert hasattr(regions, "region_element_conditions")
    assert hasattr(regions, "region_relation_conditions")
    assert not hasattr(regions, "RegionElementConditions")
    assert not hasattr(regions, "RegionRelationConditions")
    assert hasattr(regions, "is_regular_closed_region")
    assert hasattr(regions, "is_regular_open_region")


def test_redundant_inequality_does_not_create_a_singular_boundary():
    x = sp.Symbol("x", real=True)
    region = x**2 >= 0
    assert semialg.region_singular_locus(region, (x,)) is sp.false
    assert semialg.region_nonsmooth_locus(region, (x,)) is sp.false
