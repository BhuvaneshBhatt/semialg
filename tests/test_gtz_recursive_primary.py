from dataclasses import replace

import sympy as sp

from semialg import replay_certificate
from semialg.algebraic.gtz_primary import (
    gtz_primary_decomposition,
    verify_gtz_primary_decomposition_certificate,
)
from semialg.algebraic_decomposition import (
    associated_primes,
    primary_decomposition,
    verify_primary_decomposition_certificate,
)

from ._ideal_assertions import ideals_equal


def test_gtz4_nonmonomial_embedded_component():
    x, y = sp.symbols("x y")
    source = ((x + y) ** 2, y * (x + y))
    result = gtz_primary_decomposition(source, (x, y))
    assert result.complete
    assert result.irredundant
    assert len(result.components) == 2
    radicals = [component.radical for component in result.components]
    assert any(ideals_equal(radical, (x + y,), (x, y)) for radical in radicals)
    assert any(ideals_equal(radical, (x, y), (x, y)) for radical in radicals)
    assert verify_gtz_primary_decomposition_certificate(result.certificate)


def test_gtz4_same_dimension_residual_branch():
    x, y = sp.symbols("x y")
    result = gtz_primary_decomposition((x * y,), (x, y))
    assert result.complete
    assert len(result.components) == 2
    assert all(component.dimension == 1 for component in result.components)
    radicals = [component.radical for component in result.components]
    assert any(ideals_equal(radical, (x,), (x, y)) for radical in radicals)
    assert any(ideals_equal(radical, (y,), (x, y)) for radical in radicals)


def test_gtz4_localized_artinian_factorization():
    x, y = sp.symbols("x y")
    result = gtz_primary_decomposition((x * (x - 1),), (x, y))
    assert result.complete
    assert len(result.components) == 2
    radicals = [component.radical for component in result.components]
    assert any(ideals_equal(radical, (x,), (x, y)) for radical in radicals)
    assert any(ideals_equal(radical, (x - 1,), (x, y)) for radical in radicals)


def test_primary_decomposition_falls_through_to_gtz():
    x, y = sp.symbols("x y")
    source = ((x + y) ** 2, y * (x + y))
    result = primary_decomposition(source, (x, y))
    assert result.complete
    assert result.method == "gtz_recursive"
    assert verify_primary_decomposition_certificate(result.certificate)
    primes = associated_primes(source, (x, y))
    assert primes.complete
    assert len(primes.primes) == 2


def test_gtz5_unified_replay():
    x, y = sp.symbols("x y")
    result = gtz_primary_decomposition(((x + y) ** 2, y * (x + y)), (x, y))
    replay = replay_certificate(result)
    assert replay.verified is True
    public = primary_decomposition(((x + y) ** 2, y * (x + y)), (x, y))
    assert replay_certificate(public).verified is True


def test_gtz5_rejects_tampered_hull():
    x, y = sp.symbols("x y")
    result = gtz_primary_decomposition(((x + y) ** 2, y * (x + y)), (x, y))
    root = replace(result.certificate.root, hull_generators=(x,))
    forged = replace(result.certificate, root=root)
    assert not verify_gtz_primary_decomposition_certificate(forged)


def test_gtz5_rejects_tampered_final_component():
    x, y = sp.symbols("x y")
    result = gtz_primary_decomposition((x * y,), (x, y))
    first = replace(result.certificate.final_components[0], degree=2)
    forged = replace(
        result.certificate,
        final_components=(first, *result.certificate.final_components[1:]),
    )
    assert not verify_gtz_primary_decomposition_certificate(forged)


def test_gtz6_planner_and_cache_policy():
    from semialg.algebraic.gtz_primary import (
        clear_gtz_caches,
        gtz_cache_info,
        plan_gtz_primary_decomposition,
    )

    x, y, z = sp.symbols("x y z")
    plan = plan_gtz_primary_decomposition((x * y - z, x * z - y), (x, y, z))
    assert plan.variables == (x, y, z)
    assert plan.generator_count >= 1
    assert plan.dimension >= 0
    assert plan.memoize_subproblems is True
    assert plan.primitive_search == "trace_rank_then_deterministic_linear_form"

    clear_gtz_caches()
    before = gtz_cache_info()["qq"]
    gtz_primary_decomposition((x * y,), (x, y))
    middle = gtz_cache_info()["qq"]
    gtz_primary_decomposition((x * y,), (x, y))
    after = gtz_cache_info()["qq"]
    assert middle.misses > before.misses
    assert after.hits > middle.hits


def test_gtz6_result_records_execution_plan():
    x, y = sp.symbols("x y")
    result = gtz_primary_decomposition(((x + y) ** 2, y * (x + y)), (x, y))
    assert result.complete
    assert result.plan is not None
    assert result.plan.dimension == 1
    assert result.plan.memoize_subproblems
    assert verify_gtz_primary_decomposition_certificate(result.certificate)
