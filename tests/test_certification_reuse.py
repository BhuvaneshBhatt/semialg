"""Exact image and truth reuse preserves the mathematical statement."""

import sympy as sp

from semialg._optimization_certification import (
    _centered_range_problem,
    _certify_optimum_by_range,
    clear_optimization_range_cache,
    optimization_range_cache_info,
)
from semialg.context import computation_context
from semialg.optimization_results import OptimizationCertificationPolicy
from semialg.simplify.implication import (
    _shared_truth_context,
    clear_implication_minimization_stats,
    implication_minimization_stats,
)


def test_translated_image_reuse():
    x, y = sp.symbols("x y", real=True)
    policy = OptimizationCertificationPolicy("complete")
    clear_optimization_range_cache()
    original = _certify_optimum_by_range(
        x * x + y * y, sp.Eq(y * y, x * x), (x, y), kind="max", policy=policy
    )
    translated = _certify_optimum_by_range(
        (x - 3) ** 2 + (y + 2) ** 2,
        sp.Eq((y + 2) ** 2, (x - 3) ** 2),
        (x, y),
        kind="max",
        policy=policy,
    )
    assert original == translated == (sp.oo, False)
    assert optimization_range_cache_info().misses == 1
    assert optimization_range_cache_info().hits == 1


def test_translation_preserves_domain():
    x, y = sp.symbols("x y", real=True)
    objective, condition = _centered_range_problem(
        (x - 3) ** 2 + (y + 2) ** 2, sp.And(x > 3, y <= -1), (x, y)
    )
    assert objective == x * x + y * y
    for xv in (-1, 0, 1):
        for yv in (0, 1, 2):
            assert bool(condition.subs({x: xv, y: yv})) is (xv > 0 and yv <= 1)
    positive = sp.Symbol("p", positive=True)
    assert _centered_range_problem((positive - 3) ** 2, sp.true, (positive,)) == (
        (positive - 3) ** 2,
        sp.true,
    )


def test_shared_boundaries_reuse_truth():
    x, y = sp.symbols("x y", real=True)
    clear_implication_minimization_stats()
    with computation_context():
        original = _shared_truth_context(sp.And(x > 0, y > 0))
        first = original.truth_mask(x > 0)
        evaluations = implication_minimization_stats().truth_evaluations
        related = _shared_truth_context(sp.Or(x <= 0, y >= 0))
        assert related is original
        assert related.truth_mask(x > 0) == first
        assert implication_minimization_stats().truth_evaluations == evaluations
        assert implication_minimization_stats().shared_cad_builds == 1
        # Adding a boundary needs a new sign-invariant context.
        assert _shared_truth_context(sp.And(x > 0, y > 0, x + y > 0)) is not original
        assert implication_minimization_stats().shared_cad_builds == 2
    with computation_context():
        assert _shared_truth_context(sp.And(x > 0, y > 0)) is not original


def test_constant_problem_has_no_translation():
    assert _centered_range_problem(sp.S(4), sp.true, ()) == (sp.S(4), sp.true)
