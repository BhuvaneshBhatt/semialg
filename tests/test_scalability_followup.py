from __future__ import annotations

import sympy as sp


def test_shared_cad_reconstruction_merges_complementary_branches():
    from semialg.simplify.implication import (
        clear_implication_minimization_stats,
        implication_minimization_stats,
        minimize_disj_by_impl,
    )

    x, y = sp.symbols("x y", real=True)
    expr = sp.Or(sp.And(x > 0, y > 0), sp.And(x > 0, y <= 0), evaluate=False)
    clear_implication_minimization_stats()
    result = minimize_disj_by_impl(expr)
    assert result == (x > 0)
    stats = implication_minimization_stats()
    assert stats.shared_cad_builds == 1
    assert stats.cover_selected == 1
    assert stats.pairwise_fallbacks == 0


def test_linear_polytope_optimization_uses_vertex_specialization():
    from semialg import semialgebraic_minimize

    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_minimize(
        x + 2 * y,
        sp.And(x >= 0, y >= 0, x + y <= 3),
        (x, y),
        return_result=True,
    )
    assert result.value == 0
    assert result.certified
    assert result.method == "exact_linear_polytope_vertices"


def test_separable_optimization_splits_cartesian_problem():
    from semialg import semialgebraic_minimize

    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_minimize(
        (x - 1) ** 2 + (y + 2) ** 2,
        sp.And(x >= 0, x <= 4, y >= -3, y <= 1),
        (x, y),
        return_result=True,
    )
    assert result.value == 0
    assert result.method == "separable_exact_optimization"
    assert result.certified
    assert result.point == {x: 1, y: -2}


def test_public_certificate_replay_and_diagnostics():
    from semialg import replay_certificate, semialgebraic_minimize
    from semialg.certificates import result_diagnostics

    x = sp.symbols("x", real=True)
    result = semialgebraic_minimize(x, sp.And(x >= 0, x <= 2), (x,), return_result=True)
    replay = replay_certificate(result)
    diag = result_diagnostics(result)
    assert replay.verified is True
    assert diag.kind == "OptimizationResult"
    assert diag.method == result.method


def test_selective_coefficient_domain_policy():
    from semialg.cad_algorithms.polynomial_utils import (
        CoefficientDomainPolicy,
        coefficient_domain_decision,
    )

    x = sp.symbols("x")
    rational = coefficient_domain_decision(x**2 + sp.Rational(1, 3), x)
    assert rational.strategy == "natural_exact"
    algebraic = coefficient_domain_decision(x + sp.sqrt(2), x)
    assert algebraic.strategy == "small_algebraic_extension"
    guarded = coefficient_domain_decision(
        x + sp.sqrt(2),
        x,
        policy=CoefficientDomainPolicy(max_extension_degree=1),
    )
    assert guarded.strategy in {"extension_too_large", "expression_domain"}


def test_performance_probe_contracts_are_deterministic():
    from semialg.benchmarks.performance import run_core_performance_probes

    probes = {probe.name: probe for probe in run_core_performance_probes()}
    boolean = probes["shared_cad_boolean_reconstruction"]
    assert boolean.metrics["shared_cad_builds"] == 1
    assert boolean.metrics["pairwise_fallbacks"] == 0
    assert boolean.metrics["cover_selected"] == 1
    linear = probes["linear_box_optimization"]
    assert linear.metrics["method"] == "exact_linear_polytope_vertices"
    separable = probes["separable_box_optimization"]
    assert separable.metrics["method"] == "separable_exact_optimization"
