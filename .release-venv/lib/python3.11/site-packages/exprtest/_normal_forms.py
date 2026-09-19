"""Bounded shared normal forms used by multiple exact oracle stages."""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from . import _config as cfg
from ._assumptions import normalize_assumptions
from ._exact_constants import reduce_exact_trig
from ._fast import quick_reduce
from ._special import reduce_special_values
from ._transcendental import normalize_exp_log, normalize_logs


@lru_cache(maxsize=cfg.NORMAL_FORM_CACHE_SIZE)
def special_normal_form(term: sp.Expr) -> sp.Expr:
    return quick_reduce(reduce_special_values(sp.sympify(term)))


@lru_cache(maxsize=cfg.NORMAL_FORM_CACHE_SIZE)
def exp_log_normal_form(term: sp.Expr, assumptions=sp.true) -> sp.Expr:
    assumptions = normalize_assumptions(assumptions)
    return quick_reduce(normalize_exp_log(sp.sympify(term), assumptions))


@lru_cache(maxsize=cfg.NORMAL_FORM_CACHE_SIZE)
def log_normal_form(term: sp.Expr, assumptions=sp.true) -> sp.Expr:
    assumptions = normalize_assumptions(assumptions)
    return quick_reduce(normalize_logs(sp.sympify(term), assumptions))


@lru_cache(maxsize=cfg.NORMAL_FORM_CACHE_SIZE)
def trig_normal_form(term: sp.Expr) -> sp.Expr:
    return quick_reduce(reduce_exact_trig(sp.sympify(term)))


def clear_normal_form_cache() -> None:
    special_normal_form.cache_clear()
    exp_log_normal_form.cache_clear()
    log_normal_form.cache_clear()
    trig_normal_form.cache_clear()
