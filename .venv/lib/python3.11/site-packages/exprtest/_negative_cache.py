"""Bounded caches for structurally inapplicable fast-oracle stages.

Only deterministic *applicability failures* are cached here.  Mathematical
UNKNOWN verdicts are not cached as negative facts.
"""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from . import _config as cfg
from ._cost import ExactBudget, stage_allowed
from ._features import expression_features


@lru_cache(maxsize=cfg.NEGATIVE_CACHE_SIZE)
def cyclotomic_inapplicable(term: sp.Expr) -> bool:
    term = sp.sympify(term)
    if not stage_allowed(term, "cyclotomic"):
        return True
    return not expression_features(term).has_cyclotomic_shape


@lru_cache(maxsize=cfg.NEGATIVE_CACHE_SIZE)
def log_inapplicable(term: sp.Expr) -> bool:
    return not expression_features(sp.sympify(term)).has_log_exp


@lru_cache(maxsize=cfg.NEGATIVE_CACHE_SIZE)
def special_inapplicable(term: sp.Expr) -> bool:
    return not expression_features(sp.sympify(term)).has_special


@lru_cache(maxsize=cfg.NEGATIVE_CACHE_SIZE)
def log_rel_inapplicable(term: sp.Expr) -> bool:
    return not expression_features(sp.sympify(term)).has_log


@lru_cache(maxsize=cfg.NEGATIVE_CACHE_SIZE)
def tower_inapplicable(term: sp.Expr) -> bool:
    term = sp.sympify(term)
    if not stage_allowed(term, "tower"):
        return True
    features = expression_features(term)
    if features.has_tower_shape:
        return False
    return not term.has(sp.I)


@lru_cache(maxsize=cfg.NEGATIVE_CACHE_SIZE)
def pslq_inapplicable(term: sp.Expr) -> bool:
    term = sp.sympify(term)
    if not stage_allowed(term, "pslq", ExactBudget()):
        return True
    parts = term.args if term.is_Add else (term,)
    return len(parts) < 2


@lru_cache(maxsize=cfg.NEGATIVE_CACHE_SIZE)
def exp_indep_inapplicable(term: sp.Expr) -> bool:
    term = sp.sympify(term)
    features = expression_features(term)
    if not (term.is_Add and features.has_exp):
        return True
    parts = term.args
    return len(parts) < 2 or len(parts) > cfg.EXP_INDEP_MAX_TERMS


@lru_cache(maxsize=cfg.NEGATIVE_CACHE_SIZE)
def sqrt_sum_inapplicable(term: sp.Expr) -> bool:
    term = sp.sympify(term)
    features = expression_features(term)
    if not (term.is_Add and features.has_radical):
        return True
    parts = term.args
    return len(parts) < 2 or len(parts) > cfg.SQRT_SUM_MAX_TERMS


def clear_negative_cache() -> None:
    cyclotomic_inapplicable.cache_clear()
    log_inapplicable.cache_clear()
    special_inapplicable.cache_clear()
    log_rel_inapplicable.cache_clear()
    tower_inapplicable.cache_clear()
    pslq_inapplicable.cache_clear()
    exp_indep_inapplicable.cache_clear()
    sqrt_sum_inapplicable.cache_clear()
