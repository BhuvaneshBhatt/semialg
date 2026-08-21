import sympy as sp

from semialg.cad import cad_cache_stats, clear_cad_caches
from semialg.cad.projection.collins import build_collins_proj_set
from semialg.formula import parse_formula, parse_quant_form_text
from semialg.partial.qe import lazy_resolve_formula
from semialg.planner.features import extract_problem_features
from semialg.planner.heuristics import candidate_variable_orders


def test_projection_tower_cache_reuses_exact_tower():
    x, y, z = sp.symbols("x y z", real=True)
    polys = (z**2 + y * z + x, y**3 + x * y + 1, x**2 + y)
    clear_cad_caches()
    first = build_collins_proj_set(polys, (x, y, z))
    after_first = cad_cache_stats()
    second = build_collins_proj_set(polys, (x, y, z))
    after_second = cad_cache_stats()
    assert first is second
    assert after_first.projection_tower_misses == 1
    assert after_second.projection_tower_hits == 1
    assert after_second.projection_tower_misses == 1


def test_variable_order_shortlist_uses_projection_complexity_and_ec_structure():
    x, y, z = sp.symbols("x y z", real=True)
    ec = z**2 + y * z + x
    polys = (ec, y**3 + x * y + 1, x**2 + y)
    formula = parse_formula(sp.And(sp.Eq(ec, 0), y**3 + x * y + 1 >= 0, x**2 + y >= 0))
    features = extract_problem_features(formula, variables=(x, y, z))
    orders = candidate_variable_orders(features, polys, equational_constraints=(ec,), limit=6)
    assert orders
    assert orders[0].projection_poly_count is not None
    projected = [
        item.projection_poly_count for item in orders if item.projection_poly_count is not None
    ]
    assert orders[0].projection_poly_count == min(projected)
    # The EC's natural main variable z should be projected early, i.e. appear
    # at the high end of the CAD order rather than first.
    assert orders[0].order[-1] == z


def test_partial_cad_prunes_on_prefix_equational_constraint_before_descending():
    parsed = parse_quant_form_text("forall x. exists y. (x = 0) & (y^2 - 1 = 0)")
    result = lazy_resolve_formula(parsed.vars, parsed.quantifiers, parsed.matrix)
    assert result.truth_value is False
    assert result.stats.lifted_stacks == 1
    assert result.stats.visited_cells_by_level == {1: 3}
    assert result.stats.pruned_prefix_cells >= 1
    assert result.stats.ec_pruned_cells >= 1


def test_partial_cad_only_lifts_descendants_of_ec_compatible_prefix_cells():
    parsed = parse_quant_form_text("exists x. exists y. (x = 0) & (y^2 - 1 = 0)")
    result = lazy_resolve_formula(parsed.vars, parsed.quantifiers, parsed.matrix)
    assert result.truth_value is True
    # x has three cells (<0, =0, >0). Only the x=0 section needs a y stack.
    assert result.stats.lifted_stacks == 2
    assert result.stats.ec_pruned_cells >= 2
    assert result.stats.evaluated_leaf_cells == 0
