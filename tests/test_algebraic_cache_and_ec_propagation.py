import sympy as sp

from semialg import algebraic_cache_stats, clear_algebraic_caches
from semialg.algebraic.comparison import compare_samples
from semialg.algebraic.rational_univariate import compute_rational_univariate_representation
from semialg.algebraic.roots import isolate_real_roots
from semialg.algebraic.signs import sign_at_sample
from semialg.cad.bounds import AlgebraicRootFunction
from semialg.formula import parse_quant_form_text
from semialg.partial.qe import lazy_resolve_formula
from semialg.planner.features import ProblemFeatures
from semialg.planner.heuristics import candidate_variable_orders


def test_algebraic_cache_public_api_and_sign_comparison_reuse():
    clear_algebraic_caches()
    x = sp.Symbol("x", real=True)
    roots = isolate_real_roots(x**2 - 2, x)
    assert sign_at_sample(x + 3, (roots[0],)) == 1
    assert sign_at_sample(x + 3, (roots[0],)) == 1
    assert compare_samples(roots[0], roots[1]) == -1
    assert compare_samples(roots[0], roots[1]) == -1
    stats = algebraic_cache_stats()
    assert stats.sign_hits >= 1
    assert stats.comparison_hits >= 1


def test_root_function_specialization_cache_reuses_exact_result():
    clear_algebraic_caches()
    x, y = sp.symbols("x y", real=True)
    root = AlgebraicRootFunction(y**2 - x, y, 1, (x,))
    first = root.specialize({x: 4})
    before = algebraic_cache_stats()
    second = root.specialize({x: 4})
    after = algebraic_cache_stats()
    assert sp.simplify(first.as_expr() - 2) == 0
    assert sp.simplify(second.as_expr() - 2) == 0
    assert after.specialization_hits > before.specialization_hits


def test_rur_construction_cache_reuses_identical_system():
    clear_algebraic_caches()
    x, y = sp.symbols("x y", real=True)
    system = (x + y - 1, x - y)
    first = compute_rational_univariate_representation(system, (x, y))
    before = algebraic_cache_stats()
    second = compute_rational_univariate_representation(system, (x, y))
    after = algebraic_cache_stats()
    assert first == second
    assert after.rur_hits > before.rur_hits


def test_order_shortlist_reports_lifting_root_and_cell_estimates():
    x, y, z = sp.symbols("x y z", real=True)
    polys = (z**2 - x, y**2 + z - 1, x**2 - 2)
    features = ProblemFeatures(variables=(x, y, z), num_polynomials=len(polys))
    orders = candidate_variable_orders(features, polys, limit=6)
    measured = [item for item in orders if item.projection_poly_count is not None]
    assert measured
    assert all(item.estimated_lifting_roots is not None for item in measured)
    assert all(
        item.estimated_cell_count is not None and item.estimated_cell_count >= 1
        for item in measured
    )
    assert "lifting/root-count" in measured[0].reason


def test_resultant_ec_propagates_to_lower_level_and_prunes_before_z_lift():
    # z=x and z=-x are simultaneously possible only when x=0. Their resultant
    # in z is 2*x, so x=0 is a necessary derived EC before z is considered.
    parsed = parse_quant_form_text("exists x. exists z. (z - x = 0) & (z + x = 0)")
    result = lazy_resolve_formula(parsed.vars, parsed.quantifiers, parsed.matrix)
    assert result.truth_value is True
    assert result.stats.derived_ec_count >= 1
    assert result.stats.ec_section_lifts >= 1
    # Only the derived x=0 section is retained at level 1.
    assert result.stats.visited_cells_by_level[1] == 1


def test_multilevel_resultant_chain_restricts_two_prefix_levels():
    x, y, z = sp.symbols("x y z", real=True)
    # The first two ECs imply y=0 after eliminating z; y=0 together with
    # y-x=0 then implies x=0. Both lower levels can therefore section-lift.
    parsed = parse_quant_form_text(
        "exists x. exists y. exists z. (z - y = 0) & (z + y = 0) & (y - x = 0)"
    )
    result = lazy_resolve_formula(parsed.vars, parsed.quantifiers, parsed.matrix)
    assert result.truth_value is True
    assert result.stats.derived_ec_count >= 2
    assert result.stats.ec_section_lifts >= 3
    assert result.stats.visited_cells_by_level.get(1) == 1
    assert result.stats.visited_cells_by_level.get(2) == 1


def test_universal_level_is_not_unsafely_restricted_to_ec_sections():
    parsed = parse_quant_form_text("forall x. (x = 0)")
    result = lazy_resolve_formula(parsed.vars, parsed.quantifiers, parsed.matrix)
    assert result.truth_value is False
    assert result.stats.ec_section_lifts == 0
