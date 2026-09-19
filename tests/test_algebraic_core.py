import pytest
import sympy as sp

from semialg.algebraic import (
    AlgebraicRoot,
    RationalInterval,
    RationalSample,
    compare_samples,
    isolate_real_roots,
    refine_isol_intvl,
    sign_at_sample,
)
from semialg.cad_algorithms import decomp_collins_complete
from semialg.context import with_computation_context


def test_algebraic_core_01():
    with pytest.raises(TypeError):
        compare_samples(RationalSample(1), sp.Integer(2))


def test_algebraic_core_02():
    x = sp.symbols("x")
    roots = isolate_real_roots(sp.Poly(x**2 - 2, x))
    assert len(roots) == 2
    assert all(isinstance(root, AlgebraicRoot) for root in roots)
    assert compare_samples(roots[0], RationalSample(0)) < 0
    assert compare_samples(roots[1], RationalSample(0)) > 0
    assert roots[0].multiplicity == 1


def test_algebraic_core_03():
    x = sp.symbols("x")
    roots = isolate_real_roots(sp.Poly((x - 1) ** 2, x))
    assert len(roots) == 1
    assert roots[0].multiplicity == 2
    assert roots[0].interval == RationalInterval(1, 1)


def test_algebraic_core_04():
    x = sp.symbols("x")
    roots = isolate_real_roots(sp.Poly(x**2 - 2, x))
    assert sign_at_sample(sp.Poly(x**2 - 2, x), [roots[0]]) == 0
    assert sign_at_sample(sp.Poly(x, x), [roots[0]]) == -1
    assert sign_at_sample(sp.Poly(x, x), [roots[1]]) == 1


def test_algebraic_core_05():
    x = sp.symbols("x")
    root = isolate_real_roots(sp.Poly(x**2 - 2, x))[1]
    refined = refine_isol_intvl(root, steps=3)
    assert refined.polynomial == root.polynomial
    assert refined.root_index == root.root_index
    assert refined.interval.width <= root.interval.width


def test_algebraic_core_06():
    x = sp.symbols("x")
    cad = decomp_collins_complete([x**2 - 1], [x])
    assert len(cad.cells) == 5
    assert all(isinstance(cell.sample[0], (RationalSample, AlgebraicRoot)) for cell in cad.cells)
    signs = [cell.signs["x - 1"] * cell.signs["x + 1"] for cell in cad.cells]
    assert signs == [1, 0, -1, 0, 1]


def test_sector_sample_separates_overlapping_root_intervals():
    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.roots import rational_between_algebraic_reals
    from semialg.algebraic.samples import AlgebraicRoot

    x = sp.symbols("x")
    left = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1)
    right = AlgebraicRoot(sp.Poly(x**2 - 3, x), RationalInterval(1, 2), 1)
    sample = rational_between_algebraic_reals(left, right)
    assert sample**2 > 2
    assert sample**2 < 3


def test_fiber_root_refines_from_parent_certificate_without_child_expression():
    from semialg.algebraic.roots import refine_isol_intvl
    from semialg.algebraic.samples import FiberRootContext

    x, y = sp.symbols("x y", real=True)
    parent = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1)
    context = FiberRootContext(
        polynomial=sp.Poly(y**2 - x, x, y),
        variables=(x, y),
        parent_samples=(parent,),
        fiber_variable=y,
    )
    child = AlgebraicRoot(
        sp.Poly(y**2 - sp.sqrt(2), y, extension=True),
        RationalInterval(1, 2),
        1,
        root_expr=None,
        fiber_context=context,
    )

    assert child.root_expr is None
    refined = refine_isol_intvl(child, steps=3)
    assert refined.interval.width < child.interval.width
    assert refined.fiber_context == context
    assert refined.interval.left**4 < 2 < refined.interval.right**4


def test_nested_fiber_sign_uses_parent_certificates():
    from semialg.algebraic.samples import FiberRootContext
    from semialg.algebraic.signs import sign_at_sample

    x, y = sp.symbols("x y", real=True)
    parent = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1)
    context = FiberRootContext(
        polynomial=sp.Poly(y**2 - x, x, y),
        variables=(x, y),
        parent_samples=(parent,),
        fiber_variable=y,
    )
    child = AlgebraicRoot(
        sp.Poly(y**2 - sp.sqrt(2), y, extension=True),
        RationalInterval(1, 2),
        1,
        root_expr=None,
        fiber_context=context,
    )

    assert sign_at_sample(sp.Poly(y**2 - 1, x, y), (parent, child)) == 1
    assert child.root_expr is None


def test_root_identity_survives_certificate_refinement():
    from semialg.algebraic.cache import sample_certificate_key, sample_identity_key
    from semialg.algebraic.roots import refine_isol_intvl

    x = sp.symbols("x", real=True)
    root = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1)
    refined = refine_isol_intvl(root, steps=2)

    assert sample_identity_key(root) == sample_identity_key(refined)
    assert sample_certificate_key(root) != sample_certificate_key(refined)


def test_fiber_family_order_uses_root_indices():
    import semialg.algebraic.comparison as comparison
    from semialg.algebraic.samples import FiberRootContext

    x, y = sp.symbols("x y", real=True)
    parent = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1)
    context = FiberRootContext(
        polynomial=sp.Poly(y**2 - x, x, y),
        variables=(x, y),
        parent_samples=(parent,),
        fiber_variable=y,
    )
    polynomial = sp.Poly(y**2 - sp.sqrt(2), y, extension=True)
    left = AlgebraicRoot(polynomial, RationalInterval(-2, -1), 0, fiber_context=context)
    right = AlgebraicRoot(polynomial, RationalInterval(1, 2), 1, fiber_context=context)

    assert comparison.compare_samples(left, right) == -1


def test_fiber_sturm_counts_and_isolates_nested_roots():
    from semialg.algebraic.samples import FiberRootContext

    x, y = sp.symbols("x y", real=True)
    parent = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1)
    context = FiberRootContext(
        polynomial=sp.Poly(y**2 - x, x, y),
        variables=(x, y),
        parent_samples=(parent,),
        fiber_variable=y,
    )

    assert context.root_count(sp.Rational(-2), sp.Rational(2)) == 2
    assert context.root_count(sp.Rational(0), sp.Rational(2)) == 1
    isolated = context.isolate_root(1, sp.Rational(0), sp.Rational(2), steps=3)
    assert isolated.width < 1
    assert isolated.left**4 < 2 < isolated.right**4


def test_ordinary_root_refinement_reuses_sturm_certificate():
    from semialg.algebraic.cache import CACHE
    from semialg.algebraic.roots import refine_isol_intvl

    x = sp.symbols("x", real=True)
    root = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1)

    CACHE.clear()
    refined = refine_isol_intvl(root, steps=4)
    assert refined.interval.width < root.interval.width
    assert len(CACHE.root_sturm_sequences) == 1


def test_tower_escape_guard_covers_ordinary_root_inside_mixed_tower():
    from semialg.algebraic.cache import CACHE
    from semialg.algebraic.samples import FiberRootContext, tower_expression_guard

    x, y, z = sp.symbols("x y z", real=True)
    parent = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1)
    context = FiberRootContext(
        polynomial=sp.Poly(y**2 - x, x, y),
        variables=(x, y),
        parent_samples=(parent,),
        fiber_variable=y,
    )
    child = AlgebraicRoot(
        sp.Poly(y**2 - sp.sqrt(2), y, extension=True),
        RationalInterval(1, 2),
        1,
        fiber_context=context,
    )
    ordinary_outer = AlgebraicRoot(sp.Poly(z**2 - 3, z), RationalInterval(1, 2), 1)

    CACHE.stats.tower_escapes = 0
    with tower_expression_guard([parent, child, ordinary_outer]):
        ordinary_outer.as_expr()
    assert CACHE.stats.tower_escapes == 1


def test_rational_intvl_around_uses_algebraic_root_certificate():
    from semialg.algebraic.roots import rational_intvl_around

    x = sp.symbols("x", real=True)
    root = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1)

    assert rational_intvl_around(root) == root.interval


@with_computation_context
def test_fiber_family_relationship_is_cached_by_stable_parent_identity():
    from semialg.algebraic.cache import CACHE, sample_identity_key
    from semialg.algebraic.comparison import compare_samples
    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext, RationalSample

    x, y = sp.symbols("x y")
    parent = RationalSample(sp.Rational(0))
    left_context = FiberRootContext(sp.Poly(y**2 - x, x, y), (x, y), (parent,), y)
    right_context = FiberRootContext(sp.Poly((y - 1) * (y + 1), x, y), (x, y), (parent,), y)
    left = AlgebraicRoot(sp.Poly(y**2, y), RationalInterval(-1, 1), 0, fiber_context=left_context)
    right = AlgebraicRoot(
        sp.Poly(y**2 - 1, y),
        RationalInterval(sp.Rational(1, 2), sp.Rational(3, 2)),
        1,
        fiber_context=right_context,
    )
    CACHE.clear()
    # Distinct roots overlap initially; the family certificate must not infer
    # equality merely from a vanishing resultant elsewhere in the families.
    assert compare_samples(left, right) < 0
    misses = CACHE.stats.family_relationship_misses
    assert misses == 1
    # Stable parent identity, rather than interval state, owns the certificate.
    assert sample_identity_key(parent) == sample_identity_key(RationalSample(sp.Rational(0)))
    assert CACHE.stats.family_relationship_misses == misses


def test_cross_family_order_caches_sector_separator():
    from semialg.algebraic.cache import CACHE
    from semialg.algebraic.comparison import compare_samples
    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.sample_points import choose_sector_sample
    from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext, RationalSample

    x, y = sp.symbols("x y")
    parent = RationalSample(sp.Rational(0))
    left_context = FiberRootContext(sp.Poly(y**2 - 2, x, y), (x, y), (parent,), y)
    right_context = FiberRootContext(sp.Poly(y**2 - 3, x, y), (x, y), (parent,), y)
    left = AlgebraicRoot(
        sp.Poly(y**2 - 2, y), RationalInterval(1, 2), 1, fiber_context=left_context
    )
    right = AlgebraicRoot(
        sp.Poly(y**2 - 3, y), RationalInterval(1, 2), 1, fiber_context=right_context
    )
    CACHE.clear()
    assert compare_samples(left, right) < 0
    separator = choose_sector_sample(left, right)
    assert left.interval.left < separator.value < right.interval.right
    assert CACHE.stats.comparison_refinements == 0


def test_cross_family_common_root_uses_subresultant_certificate():
    from semialg.algebraic.cache import CACHE
    from semialg.algebraic.comparison import compare_samples
    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext, RationalSample

    x, y = sp.symbols("x y")
    parent = RationalSample(sp.Rational(0))
    left_context = FiberRootContext(sp.Poly(y**2 - x, x, y), (x, y), (parent,), y)
    right_context = FiberRootContext(sp.Poly(y * (y - 1), x, y), (x, y), (parent,), y)
    left = AlgebraicRoot(sp.Poly(y**2, y), RationalInterval(-1, 1), 0, fiber_context=left_context)
    right = AlgebraicRoot(
        sp.Poly(y * (y - 1), y),
        RationalInterval(-1, sp.Rational(1, 2)),
        0,
        fiber_context=right_context,
    )
    CACHE.clear()
    assert compare_samples(left, right) == 0
    assert CACHE.stats.comparison_refinements == 0


def test_sparse_tower_remainder_matches_exact_triangular_reduction():
    from semialg.algebraic.signs import _monic_sparse_remainder, _sparse_qq_expr, _sparse_qq_poly

    x, y, z = sp.symbols("x y z")
    gens = (x, y, z)
    target_expr = z**4 + x * z**3 + y * z**2 + x * y * z + 3
    defining_expr = 2 * z**2 + x * z + y
    target = _sparse_qq_poly(target_expr, gens)
    defining = _sparse_qq_poly(defining_expr, gens)
    assert target is not None and defining is not None

    sparse_remainder = _monic_sparse_remainder(target, defining, 2)
    assert sparse_remainder is not None
    expected = sp.Poly(target_expr, z, domain=sp.QQ.poly_ring(x, y)).rem(
        sp.Poly(defining_expr, z, domain=sp.QQ.poly_ring(x, y))
    )
    assert sp.expand(_sparse_qq_expr(sparse_remainder, gens) - expected.as_expr()) == 0


@with_computation_context
def test_fiber_root_index_partition_reuses_certified_splits():
    from semialg.algebraic.cache import CACHE
    from semialg.algebraic.samples import FiberRootContext, RationalSample

    x, y = sp.symbols("x y", real=True)
    context = FiberRootContext(
        polynomial=sp.Poly(y**3 - 2 * y, x, y),
        variables=(x, y),
        parent_samples=(RationalSample(0),),
        fiber_variable=y,
    )

    CACHE.clear()
    first = context.isolate_root(2, sp.Rational(-2), sp.Rational(2), steps=2)
    split_count = len(CACHE.root_intvl_splits)
    endpoint_count = len(CACHE.sturm_endpoint_certs)
    repeated = context.isolate_root(2, sp.Rational(-2), sp.Rational(2), steps=0)

    assert repeated == first
    assert context.best_root_intvl(2) == first
    assert len(CACHE.root_intvl_splits) == split_count
    assert len(CACHE.sturm_endpoint_certs) == endpoint_count


def test_sign_interval_certificate_is_shared_across_stale_root_copies():
    from semialg.algebraic.cache import CACHE, clear_algebraic_caches
    from semialg.algebraic.roots import publish_sign_intvl, refresh_sign_root

    x = sp.symbols("x", real=True)
    clear_algebraic_caches()
    broad = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1)
    narrow = AlgebraicRoot(
        sp.Poly(x**2 - 2, x), RationalInterval(sp.Rational(7, 5), sp.Rational(3, 2)), 1
    )

    publish_sign_intvl(broad)
    publish_sign_intvl(narrow)
    refreshed = refresh_sign_root(broad)

    assert refreshed.interval == narrow.interval
    assert len(CACHE.sign_intvl_certs) == 1


@with_computation_context
def test_fiber_root_rank_cert_avoids_recounting_known_rank_boundaries():
    from semialg.algebraic.cache import CACHE
    from semialg.algebraic.samples import FiberRootContext, RationalSample

    x, y = sp.symbols("x y", real=True)
    context = FiberRootContext(
        polynomial=sp.Poly((y + 1) * y * (y - 1), x, y),
        variables=(x, y),
        parent_samples=(RationalSample(0),),
        fiber_variable=y,
    )

    CACHE.clear()
    roots = context.isolate_all_roots()
    cert = context.root_rank_cert(2)
    assert cert is not None
    assert cert.root_index == 2
    assert cert.node_root_count == 1
    endpoint_count = len(CACHE.sturm_endpoint_certs)

    repeated = context.isolate_root(2, sp.Rational(-2), sp.Rational(2), steps=0)
    assert repeated == roots[2]
    assert context.root_on_or_below(2, cert.interval.left - 1) is False
    assert context.root_on_or_below(2, cert.interval.right + 1) is True
    assert len(CACHE.sturm_endpoint_certs) == endpoint_count


@with_computation_context
def test_ranked_partition_reuses_existing_rational_boundary_without_sturm_recount():
    from semialg.algebraic.cache import CACHE
    from semialg.algebraic.samples import FiberRootContext, RationalSample

    x, y = sp.symbols("x y", real=True)
    context = FiberRootContext(
        polynomial=sp.Poly(y**3 - 2 * y, x, y),
        variables=(x, y),
        parent_samples=(RationalSample(0),),
        fiber_variable=y,
    )

    CACHE.clear()
    context.isolate_all_roots()
    assert len(CACHE.ranked_partition_nodes) > 0
    assert len(CACHE.boundary_root_ranks) > 0
    endpoint_count = len(CACHE.sturm_endpoint_certs)

    assert context.roots_left_of(sp.Rational(0)) == 2
    assert context.root_on_or_below(0, sp.Rational(0)) is True
    assert context.root_on_or_below(2, sp.Rational(0)) is False
    assert len(CACHE.sturm_endpoint_certs) == endpoint_count


def test_descartes_variation_and_split_nodes_are_persistent():
    from semialg.algebraic.cache import CACHE
    from semialg.algebraic.roots import _admissible_descartes_split, _descartes_variations_interval

    CACHE.clear()
    x = sp.symbols("x")
    poly = sp.Poly(x**5 - 3 * x + 1, x, domain=sp.QQ)
    left, right = sp.Rational(-3), sp.Rational(3)
    variations = _descartes_variations_interval(poly, left, right)
    split = _admissible_descartes_split(poly, left, right, variations)
    variation_count = len(CACHE.descartes_variations)
    node_count = len(CACHE.descartes_nodes)
    assert _descartes_variations_interval(poly, left, right) == variations
    assert _admissible_descartes_split(poly, left, right, variations) == split
    assert len(CACHE.descartes_variations) == variation_count
    assert len(CACHE.descartes_nodes) == node_count


@with_computation_context
def test_root_count_promotes_sturm_endpoints_to_ranked_boundary_certificates():
    from semialg.algebraic.cache import CACHE
    from semialg.algebraic.samples import FiberRootContext, RankedBoundaryCert, RationalSample

    x, y = sp.symbols("x y", real=True)
    context = FiberRootContext(
        polynomial=sp.Poly(y**3 - 2 * y, x, y),
        variables=(x, y),
        parent_samples=(RationalSample(0),),
        fiber_variable=y,
    )
    CACHE.clear()
    assert context.root_count(sp.Rational(-3, 2), sp.Rational(3, 2)) == 3
    left = CACHE.ranked_boundary_certs.get(context._boundary_cert_key(sp.Rational(-3, 2)))
    right = CACHE.ranked_boundary_certs.get(context._boundary_cert_key(sp.Rational(3, 2)))
    assert isinstance(left, RankedBoundaryCert) and left.sturm is not None
    assert isinstance(right, RankedBoundaryCert) and right.sturm is not None
    endpoint_count = len(CACHE.sturm_endpoint_certs)
    assert context.roots_left_of(sp.Rational(-3, 2)) == left.root_rank
    assert context.roots_left_of(sp.Rational(3, 2)) == right.root_rank
    assert len(CACHE.sturm_endpoint_certs) == endpoint_count


def test_fiber_zero_uses_nonzero_constant_subresultant_as_coprime_certificate():
    from semialg.algebraic.cache import CACHE
    from semialg.algebraic.comparison import (
        _fiber_common_factor_cache_key,
        _fiber_poly_zero_at_root,
    )
    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext, RationalSample

    x, y = sp.symbols("x y")
    parent = RationalSample(sp.Rational(0))
    defining_context = FiberRootContext(sp.Poly(y**2 - 2, x, y), (x, y), (parent,), y)
    root = AlgebraicRoot(
        sp.Poly(y**2 - 2, y), RationalInterval(1, 2), 1, fiber_context=defining_context
    )
    target = sp.Poly(y - 1, x, y)
    target_context = FiberRootContext(target, (x, y), (parent,), y)

    CACHE.clear()
    assert not _fiber_poly_zero_at_root(target, root)
    cache_key = _fiber_common_factor_cache_key(defining_context, target_context)
    assert CACHE.fiber_common_gcds.get(cache_key) is False


def test_fiber_context_detects_section_nullification_before_root_isolation():
    """A generic fiber root must disappear when every coefficient vanishes on the section."""
    import sympy as sp

    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext, RationalSample

    y, x, t = sp.symbols("y x t")
    x_poly = sp.Poly(x**2 + 14 * x + 44, x, domain=sp.QQ)
    x_root = AlgebraicRoot(
        x_poly,
        RationalInterval(sp.Rational(-10), sp.Rational(-9)),
        0,
        root_expr=-7 - sp.sqrt(5),
    )
    context = FiberRootContext(
        sp.Poly(t * ((x - y) ** 2 - 5), y, x, t, domain=sp.QQ),
        (y, x, t),
        (RationalSample(-7), x_root),
        t,
    )

    assert context.is_nullified_at_parent()


def test_local_tower_root_expr_rejects_nullified_fiber_root():
    import sympy as sp

    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.samples import (
        AlgebraicRoot,
        FiberRootContext,
        RationalSample,
        local_tower_root_expr,
    )

    y, x, t = sp.symbols("y x t")
    x_poly = sp.Poly(x**2 + 14 * x + 44, x, domain=sp.QQ)
    x_root = AlgebraicRoot(
        x_poly,
        RationalInterval(sp.Rational(-10), sp.Rational(-9)),
        0,
        root_expr=-7 - sp.sqrt(5),
    )
    context = FiberRootContext(
        sp.Poly(t * ((x - y) ** 2 - 5), y, x, t, domain=sp.QQ),
        (y, x, t),
        (RationalSample(-7), x_root),
        t,
    )
    stale_generic_root = AlgebraicRoot(
        sp.Poly(t, t, domain=sp.QQ), RationalInterval(0, 0), 0, fiber_context=context
    )

    assert local_tower_root_expr(stale_generic_root) is None


def test_sign_eliminates_certified_point_coordinate_before_tower_interval_arithmetic():
    import sympy as sp

    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext
    from semialg.algebraic.signs import sign_at_sample

    y, x = sp.symbols("y x")
    y_root = AlgebraicRoot(
        sp.Poly(4 * y**2 - 3, y, domain=sp.QQ),
        RationalInterval(sp.Rational(-1), sp.Rational(-3, 4)),
        0,
    )
    context = FiberRootContext(
        sp.Poly(y**2 + x**2 - 1, y, x, domain=sp.QQ),
        (y, x),
        (y_root,),
        x,
    )
    x_root = AlgebraicRoot(
        sp.Poly(x**2 + y**2 - 1, x, domain=sp.QQ.frac_field(y)),
        RationalInterval(sp.Rational(-1, 2), sp.Rational(-1, 2)),
        0,
        fiber_context=context,
    )

    assert sign_at_sample(sp.Poly(-2 * x * y - y, y, x, domain=sp.QQ), (y_root, x_root)) == 0


def test_multilevel_quotient_sturm_avoids_global_zero_divisor_inverse():
    """The selected tower is a field even when its QQ quotient has conjugate zero divisors."""
    import sympy as sp

    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext

    x, y, t = sp.symbols("x y t")
    x_root = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1)
    y_context = FiberRootContext(sp.Poly(y**2 - 2, x, y), (x, y), (x_root,), y)
    y_root = AlgebraicRoot(sp.Poly(y**2 - 2, y), RationalInterval(1, 2), 1, fiber_context=y_context)
    # QQ[x,y]/(x^2-2,y^2-2) has zero divisors: x+y vanishes on a conjugate
    # component.  It is nevertheless nonzero at the selected (+sqrt(2),+sqrt(2))
    # section, so a specialization-safe PRS must not try to invert it globally.
    context = FiberRootContext(
        sp.Poly((x + y) * t**2 + t + 1, x, y, t),
        (x, y, t),
        (x_root, y_root),
        t,
    )

    assert context._tower_quotient_sturm_sequence is not None
    assert context.root_count(sp.Rational(-10), sp.Rational(10)) == 0


def test_multilevel_quotient_sturm_eliminates_exact_parent_coordinates():
    import sympy as sp

    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext, RationalSample

    y, x, c, t = sp.symbols("y x c t")
    y_sample = RationalSample(sp.Rational(-1))
    x_context = FiberRootContext(sp.Poly(y**2 + x**2 - 3, y, x), (y, x), (y_sample,), x)
    x_root = AlgebraicRoot(sp.Poly(x**2 - 2, x), RationalInterval(1, 2), 1, fiber_context=x_context)
    c_context = FiberRootContext(sp.Poly(c**2 - x**2, y, x, c), (y, x, c), (y_sample, x_root), c)
    c_root = AlgebraicRoot(
        sp.Poly(c**2 - x**2, c, domain=sp.QQ.frac_field(x)),
        RationalInterval(1, 2),
        1,
        fiber_context=c_context,
    )
    context = FiberRootContext(
        sp.Poly((x + c) * t**2 + t + 1, y, x, c, t),
        (y, x, c, t),
        (y_sample, x_root, c_root),
        t,
    )

    assert context._tower_quotient_sturm_sequence is not None
    assert context.root_count(sp.Rational(-10), sp.Rational(10)) == 0


def test_projected_sign_reuses_selected_fiber_zero_certificate_with_omitted_parent():
    """A sign polynomial may omit an algebraic parent carried by the final fiber root."""
    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext, RationalSample
    from semialg.algebraic.signs import sign_at_sample

    x, y, c, b = sp.symbols("x y c b")
    x_sample = RationalSample(sp.Rational(-1, 2))
    y_context = FiberRootContext(sp.Poly(x**2 + y**2 - 1, x, y), (x, y), (x_sample,), y)
    y_root = AlgebraicRoot(
        sp.Poly(y**2 - sp.Rational(3, 4), y),
        RationalInterval(sp.Rational(-7, 8), sp.Rational(-3, 4)),
        0,
        fiber_context=y_context,
    )
    c_context = FiberRootContext(
        sp.Poly(x**2 - 2 * x * c + c**2 - 1, x, y, c),
        (x, y, c),
        (x_sample, y_root),
        c,
    )
    c_root = AlgebraicRoot(
        sp.Poly(c**2 + c - sp.Rational(3, 4), c),
        RationalInterval(sp.Rational(0), sp.Rational(1)),
        1,
        fiber_context=c_context,
    )
    b_context = FiberRootContext(
        sp.Poly(x * y**2 * b - x * b**3 - y**2 * c * b + c * b**3, x, y, c, b),
        (x, y, c, b),
        (x_sample, y_root, c_root),
        b,
    )
    b_root = AlgebraicRoot(
        sp.Poly((c - x) * b**3 + (-c * y**2 + x * y**2) * b, b, domain=sp.QQ.frac_field(x, y, c)),
        RationalInterval(sp.Rational(-1), sp.Rational(-1, 2)),
        0,
        fiber_context=b_context,
    )

    assert sign_at_sample(sp.Poly(-(b**2) + y**2, y, b), (y_root, b_root)) == 0


def test_fiber_context_nullification_accepts_algebraic_extension_coefficients():
    """Section nullification must not force algebraic coefficients into QQ(parent)."""
    import sympy as sp

    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext

    z, t = sp.symbols("z t")
    parent = AlgebraicRoot(
        sp.Poly(z**2 - 3, z, domain=sp.QQ),
        RationalInterval(sp.Rational(1), sp.Rational(2)),
        1,
        root_expr=sp.sqrt(3),
    )
    nullified = FiberRootContext(
        sp.Poly(t * (z - sp.sqrt(3)), z, t, extension=True),
        (z, t),
        (parent,),
        t,
    )
    nonnullified = FiberRootContext(
        sp.Poly(t * (z + sp.sqrt(3)), z, t, extension=True),
        (z, t),
        (parent,),
        t,
    )

    assert nullified.is_nullified_at_parent()
    assert not nonnullified.is_nullified_at_parent()
