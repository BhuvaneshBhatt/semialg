"""High-level exact algebraic zero-testing operations."""

from __future__ import annotations

import sympy as sp

from ._algebraic_bounds import (
    AlgebraicGap,
    algebraic_gap_bound,
    algebraic_gap_test,
)
from ._algebraic_model import (
    AlgebraicExpressionInfo,
    AlgebraicModel,
    algebraic_relation_reduction_test,
    build_algebraic_model,
    classify_algebraic_expression,
    common_number_field_test,
    extract_algebraic_generators,
)
from ._algebraic_normalize import normalize_radicals, simplify_root_sets
from ._algebraic_tower import (
    AlgebraicTower,
    TowerStep,
    build_algebraic_tower,
    tower_algebraic_test,
)
from ._cost import within_budget
from ._result import Verdict, ZeroClassification

__all__ = [
    "AlgebraicExpressionInfo",
    "AlgebraicGap",
    "AlgebraicModel",
    "AlgebraicTower",
    "TowerStep",
    "algebraic_gap_bound",
    "algebraic_gap_test",
    "algebraic_relation_reduction_test",
    "build_algebraic_model",
    "build_algebraic_tower",
    "classify_algebraic_expression",
    "common_number_field_test",
    "exact_algebraic_number_test",
    "extract_algebraic_generators",
    "normalize_radicals",
    "simplify_root_sets",
    "tower_algebraic_test",
]


def exact_algebraic_number_test(expr: sp.Expr) -> ZeroClassification:
    """Proof-oriented exact test for closed algebraic expressions."""
    if not within_budget(sp.sympify(expr)):
        return ZeroClassification(
            Verdict.UNKNOWN, "algebraic", detail="exact-method budget exceeded"
        )

    expr = sp.sympify(expr)
    info = classify_algebraic_expression(expr)
    if not info.is_algebraic:
        return ZeroClassification(Verdict.UNKNOWN, "algebraic", detail=info.detail)

    reduced = algebraic_relation_reduction_test(expr)
    if reduced.verdict is not Verdict.UNKNOWN:
        return reduced

    from ._exact_constants import canonical_algebraic, reduce_exact_trig

    prepared = reduce_exact_trig(expr)
    prepared = normalize_radicals(simplify_root_sets(prepared))
    compact = canonical_algebraic(prepared)
    if compact == 0:
        return ZeroClassification(
            Verdict.ZERO_PROVEN,
            "algebraic-canonical",
            detail="exact algebraic canonicalization reduced expression to zero",
            evidence="exact-algebraic-canonicalization",
        )
    if compact != prepared and classify_algebraic_expression(compact).is_algebraic:
        prepared = compact
    if prepared != expr:
        reduced = algebraic_relation_reduction_test(prepared)
        if reduced.verdict is not Verdict.UNKNOWN:
            return ZeroClassification(
                reduced.verdict,
                "algebraic-preprocess",
                detail=f"exact root/radical normalization; {reduced.detail}",
                evidence=reduced.evidence,
            )

    # Primitive-element conversion is usually cheaper for one or two flat
    # generators; sparse tower arithmetic wins more often for nested/multiple
    # radicals.  Choose rather than always paying both in the same order.
    generators = extract_algebraic_generators(prepared)
    common_first = len(generators) <= 2 and not any(
        g.is_Pow and g.base.has(sp.Pow) for g in generators
    )
    methods = (
        (common_number_field_test, tower_algebraic_test)
        if common_first
        else (tower_algebraic_test, common_number_field_test)
    )
    for method in methods:
        result = method(prepared)
        if result.verdict is not Verdict.UNKNOWN:
            return result

    from ._exact_constants import algebraic_sign

    sign = algebraic_sign(prepared)
    if sign is not None:
        return ZeroClassification(
            Verdict.ZERO_PROVEN if sign == 0 else Verdict.NONZERO_PROVEN,
            "algebraic-order",
            detail=(
                "exact real-algebraic comparison identified zero"
                if sign == 0
                else "exact real-algebraic comparison established a strict sign"
            ),
            evidence="exact-real-algebraic-order",
        )

    gap = algebraic_gap_test(prepared)
    if gap.verdict is not Verdict.UNKNOWN:
        return gap
    return ZeroClassification(
        Verdict.UNKNOWN, "algebraic", detail="exact algebraic methods were inconclusive"
    )
