"""Bounded cache for proof-grade zero classifications."""

from __future__ import annotations

from collections import OrderedDict

import sympy as sp

from . import _config as cfg
from ._result import ZeroClassification

_RESULT_CACHE: OrderedDict = OrderedDict()


def make_cache_key(expr: sp.Expr, assumptions: sp.Expr, sz_target_error: float):
    """Build a cache key from deterministic options affecting classification."""
    return (expr, assumptions, float(sz_target_error))


def cache_get(key):
    value = _RESULT_CACHE.get(key)
    if value is not None:
        _RESULT_CACHE.move_to_end(key)
    return value


def cache_set(key, value: ZeroClassification):
    """Cache only proof-grade results; never freeze heuristic/random outcomes."""
    if not value.proven:
        return
    _RESULT_CACHE[key] = value
    _RESULT_CACHE.move_to_end(key)
    while len(_RESULT_CACHE) > cfg.RESULT_CACHE_SIZE:
        _RESULT_CACHE.popitem(last=False)


def cache_clear() -> None:
    """Clear result, structural, normalization, and exact-object caches."""
    from ._exp_independence import clear_exp_indep_cache
    from ._features import clear_feature_cache
    from ._memo import clear_exact_cache
    from ._negative_cache import clear_negative_cache
    from ._nonzero import clear_nonzero_cache
    from ._normal_forms import clear_normal_form_cache
    from ._rewrite import clear_rewrite_cache
    from ._sqrt_sum import clear_sqrt_sum_cache

    _RESULT_CACHE.clear()
    clear_exact_cache()
    clear_exp_indep_cache()
    clear_sqrt_sum_cache()
    clear_negative_cache()
    clear_feature_cache()
    clear_nonzero_cache()
    clear_normal_form_cache()
    clear_rewrite_cache()
