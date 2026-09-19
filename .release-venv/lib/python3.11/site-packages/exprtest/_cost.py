"""Structural admission budgets for exact zero-testing methods."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from . import _config as cfg
from ._features import ExprFeatures, expression_features


@dataclass(frozen=True)
class ExactBudget:
    """Structural limits used before potentially expensive exact operations."""

    max_ops: int = cfg.EXACT_MAX_OPS
    max_nodes: int = cfg.EXACT_MAX_NODES
    max_depth: int = cfg.EXACT_MAX_DEPTH
    max_gens: int = cfg.EXACT_MAX_GENERATORS
    max_degree: int = cfg.EXACT_MAX_DEGREE_PRODUCT
    max_cyclo: int = cfg.MAX_CYCLOTOMIC_ORDER
    max_logs: int = cfg.EXACT_MAX_LOG_TERMS
    max_pow_bits: int = cfg.EXACT_MAX_POW_BITS
    max_int_bits: int = cfg.EXACT_MAX_INT_BITS
    max_poly_terms: int = cfg.EXACT_MAX_POLY_TERMS
    max_result_vars: int = cfg.EXACT_MAX_RESULTANT_VARS
    max_pslq_terms: int = cfg.EXACT_MAX_PSLQ_TERMS
    pslq_digits: int = cfg.EXACT_PSLQ_PREC_DIGITS
    pslq_coeff: int = cfg.EXACT_PSLQ_MAX_COEFF


def within_budget(
    term: sp.Expr,
    budget: ExactBudget | None = None,
    *,
    features: ExprFeatures | None = None,
) -> bool:
    """Return whether ``term`` fits the common exact-method budget."""
    active = budget or ExactBudget()
    feats = features or expression_features(sp.sympify(term))
    return (
        feats.nodes <= active.max_nodes
        and feats.depth <= active.max_depth
        and feats.generators <= active.max_gens
        and feats.ops <= active.max_ops
        and feats.max_pow_bits <= active.max_pow_bits
        and feats.max_int_bits <= active.max_int_bits
    )


def stage_allowed(
    term: sp.Expr,
    stage: str,
    budget: ExactBudget | None = None,
    *,
    features: ExprFeatures | None = None,
) -> bool:
    """Apply a method-specific gate using one structural fingerprint."""
    active = budget or ExactBudget()
    feats = features or expression_features(sp.sympify(term))
    if not within_budget(term, active, features=feats):
        return False
    if stage == "resultant":
        return feats.generators <= active.max_result_vars and feats.ops <= active.max_ops // 2
    if stage == "tower":
        return (
            feats.generators <= min(active.max_gens, cfg.TOWER_MAX_GENS)
            and feats.ops <= min(active.max_ops, cfg.EXACT_TOWER_MAX_OPS)
        )
    if stage == "pslq":
        parts = len(term.args) if term.is_Add else 1
        return parts <= active.max_pslq_terms and feats.ops <= active.max_ops // 2
    if stage == "minpoly":
        return feats.generators <= active.max_gens and feats.ops <= active.max_ops
    if stage == "cyclotomic":
        return feats.ops <= active.max_ops // 2
    if stage == "factor-rational":
        return feats.max_int_bits <= cfg.EXACT_MAX_FACTOR_BITS and feats.logs <= active.max_logs
    if stage == "special":
        return feats.depth <= cfg.EXACT_MAX_SPECIAL_DEPTH
    return True
