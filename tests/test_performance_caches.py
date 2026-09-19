import sympy as sp

from semialg.cad_algorithms import clear_cad_caches
from semialg.cad_algorithms.decomposition import decomp_collins_complete
from semialg.cad_algorithms.polynomial_utils import (
    _resultant_expression_cached,
    exact_univariate_poly,
    normalize_poly,
    polynomial_key,
)
from semialg.planner.features import ProblemFeatures
from semialg.planner.heuristics import (
    _count_distinct_real_roots,
    _probe_univariate_poly,
    candidate_variable_orders,
    clear_variable_order_cache,
    variable_order_cache_info,
)
from semialg.simplify.implication import (
    clear_implication_minimization_stats,
    implication_minimization_stats,
    minimize_disj_by_impl,
)


def test_polynomial_key_memoizes_readable_key() -> None:
    x, y = sp.symbols("x y", real=True)
    poly = sp.Poly(x**2 + x * y + 1, x, y)
    polynomial_key.cache_clear()
    first = polynomial_key(poly)
    before = polynomial_key.cache_info()
    second = polynomial_key(poly)
    after = polynomial_key.cache_info()
    assert first == second
    assert after.hits == before.hits + 1


def test_projection_operation_caches_reuse_exact_work() -> None:
    x, y = sp.symbols("x y", real=True)
    left = sp.Poly(y**3 + x * y + 1, x, y)
    right = sp.Poly(y**2 + x + 1, x, y)

    normalize_poly.cache_clear()
    first_normalized = normalize_poly(left)
    normalize_before = normalize_poly.cache_info()
    second_normalized = normalize_poly(left)
    normalize_after = normalize_poly.cache_info()
    assert first_normalized is second_normalized
    assert normalize_after.hits == normalize_before.hits + 1

    # Test the memoized fallback at the function that owns that cache. FLINT
    # acceleration is a separate path and intentionally bypasses this cache.
    _resultant_expression_cached.cache_clear()
    first_resultant = _resultant_expression_cached(left.as_expr(), right.as_expr(), y)
    resultant_before = _resultant_expression_cached.cache_info()
    second_resultant = _resultant_expression_cached(left.as_expr(), right.as_expr(), y)
    resultant_after = _resultant_expression_cached.cache_info()
    assert sp.expand(first_resultant - second_resultant) == 0
    assert resultant_after.hits == resultant_before.hits + 1


def test_complete_cad_records_sign_invariance_diagnostics() -> None:
    x, y = sp.symbols("x y", real=True)
    clear_cad_caches()
    cad = decomp_collins_complete((x**2 + y**2 - 1, x + y), (x, y))
    assert cad.diagnostics is not None
    assert cad.diagnostics.invariant_failures == ()


def test_probe_poly_prefers_natural_exact_domains() -> None:
    x = sp.Symbol("x", real=True)
    assert _probe_univariate_poly(x**2 + 2 * x + 1, x).domain == sp.ZZ
    assert _probe_univariate_poly(x**2 + sp.Rational(1, 2), x).domain == sp.QQ
    algebraic = _probe_univariate_poly(x**2 + sp.sqrt(2), x)
    assert algebraic is not None
    assert algebraic.domain != sp.EX


def test_candidate_variable_orders_caches_complete_result() -> None:
    x, y = sp.symbols("x y", real=True)
    polys = [y**2 - x, x**2 + y]
    features = ProblemFeatures(variables=(x, y), num_polynomials=2)
    clear_variable_order_cache()
    first = candidate_variable_orders(features, polys, limit=4)
    before = variable_order_cache_info()
    second = candidate_variable_orders(features, list(polys), limit=4)
    after = variable_order_cache_info()
    assert first is second
    assert after.hits == before.hits + 1


def test_implication_minimization_profiling_counters() -> None:
    x, y = sp.symbols("x y", real=True)
    expr = sp.Or(sp.And(x > 1, y > 0), sp.And(x >= 0, y > 0))
    clear_implication_minimization_stats()
    result = minimize_disj_by_impl(expr)
    stats = implication_minimization_stats()
    assert result == sp.And(x >= 0, y > 0)
    assert stats.minimize_calls == 1
    assert stats.input_branches >= 2
    assert stats.output_branches == 1
    assert stats.conjunction_checks > 0
    assert stats.disjunction_checks > 0
    assert stats.cad_requests > 0
    assert stats.shared_cad_builds == 1
    assert stats.truth_evaluations > 0
    assert stats.pairwise_fallbacks == 0
    assert stats.cad_requests == 1


def test_shared_cad_implication_minimizer_uses_one_decomposition() -> None:
    x = sp.Symbol("x", real=True)
    expr = sp.Or(x > 2, x > 1, x >= 0)
    clear_implication_minimization_stats()
    result = minimize_disj_by_impl(expr)
    stats = implication_minimization_stats()
    assert result == (x >= 0)
    assert stats.shared_cad_builds == 1
    assert stats.cad_requests == 1
    assert stats.removed_branches == 2


def test_distinct_root_count_avoids_constructing_duplicate_roots() -> None:
    x = sp.Symbol("x", real=True)
    polys = (sp.Poly((x - 1) ** 2 * (x + 2), x), sp.Poly((x - 1) * (x - 3), x))
    assert _count_distinct_real_roots(polys, x) == 3


def test_exact_univariate_poly_preserves_ground_and_algebraic_domains() -> None:
    x = sp.Symbol("x", real=True)
    assert exact_univariate_poly(x**2 - 2, x).domain == sp.ZZ
    assert exact_univariate_poly(x**2 - sp.Rational(1, 2), x).domain == sp.QQ
    algebraic = exact_univariate_poly(x**2 - sp.sqrt(2), x)
    assert algebraic.domain != sp.EX
    assert getattr(algebraic.domain, "is_Exact", False)
    # Inexact coefficients use the conservative EX domain.
    assert exact_univariate_poly(x**2 - sp.Float("0.5"), x).domain == sp.EX


def test_adaptive_pilot_can_skip_uncompetitive_runner_up() -> None:
    x, y, z = sp.symbols("x y z", real=True)
    polys = (z**6 + x**2, y + 1, x + 1)
    features = ProblemFeatures(variables=(x, y, z), num_polynomials=len(polys))
    clear_variable_order_cache()
    scores = candidate_variable_orders(features, polys, limit=6)
    assert scores
    assert sum(score.pilot_lifting_roots is not None for score in scores) == 1


def test_optimization_range_certification_cache_reuses_complete_image() -> None:
    from semialg._optimization_certification import _cached_complete_range_certification
    from semialg.optimization import (
        clear_optimization_range_cache,
        optimization_range_cache_info,
    )

    x = sp.Symbol("x", real=True)
    clear_optimization_range_cache()
    objective = x**2
    condition = sp.And(x >= -2, x <= 3)
    first = _cached_complete_range_certification(objective, condition, (x,))
    before = optimization_range_cache_info()
    second = _cached_complete_range_certification(objective, condition, (x,))
    after = optimization_range_cache_info()
    assert first is not None
    assert first is second
    assert after.hits == before.hits + 1
