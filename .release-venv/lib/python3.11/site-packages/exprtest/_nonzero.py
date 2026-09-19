"""Very cheap structural nonzero proofs for the fast oracle."""

from __future__ import annotations

from collections.abc import MutableMapping
from functools import lru_cache
from typing import TYPE_CHECKING

import sympy as sp

from . import _config as cfg
from ._assumptions import (
    ask_property,
    assumption_facts,
    is_nonzero,
    is_positive,
    normalize_assumptions,
)
from ._defined import quick_defined, quick_defined_nonzero
from ._domains import NumberKind, number_kind
from ._errors import RECOVERABLE_SYMPY_ERRORS
from ._exp_independence import exponential_independence_nonzero
from ._fast import literal_number_zero
from ._features import expression_features
from ._function_facts import PROPERTY_FUNCS, function_nonzero
from ._negative_cache import exp_indep_inapplicable

if TYPE_CHECKING:
    from ._call_memo import OracleMemo


def _ask_nonzero(term: sp.Expr, assumptions) -> bool | None:
    try:
        facts = assumption_facts(assumptions)
        if facts:
            with sp.assumptions.assuming(*facts):
                value = sp.ask(sp.Q.nonzero(term))
        else:
            value = sp.ask(sp.Q.nonzero(term))
        return None if value is None else bool(value)
    except RECOVERABLE_SYMPY_ERRORS:
        return None


def _gamma_nonzero(term: sp.Expr) -> bool | None:
    if term.func is not sp.gamma or len(term.args) != 1:
        return None
    arg = term.args[0]
    # Gamma has no zeros. At poles the expression is not a finite nonzero
    # number, so only use the theorem when a pole can be excluded cheaply.
    if arg.is_integer is True and arg.is_nonpositive is True:
        return None
    if arg.is_finite is True and (arg.is_positive is True or arg.is_integer is False):
        return True
    return None


@lru_cache(maxsize=cfg.EXACT_NONZERO_CACHE_SIZE)
def _closed_nonzero(term: sp.Expr) -> bool | None:
    """Cached closed-expression portion of :func:`quick_nonzero`."""
    term = sp.sympify(term)
    features = expression_features(term)
    if features.has_float or features.has_nonfinite:
        return None
    if features.nodes <= 12 and term.is_zero is True:
        return False
    literal_zero = literal_number_zero(term)
    if literal_zero is not None:
        return not literal_zero
    if term.is_Float:
        return None
    if term in (sp.pi, sp.E, sp.I, sp.GoldenRatio):
        return True
    if term.func is sp.exp and len(term.args) == 1:
        arg = term.args[0]
        if arg.is_finite is True:
            return True
    gamma = _gamma_nonzero(term)
    if gamma is not None:
        return gamma
    if term.is_Mul:
        values = tuple(_closed_nonzero(arg) for arg in term.args)
        if all(value is True for value in values):
            return True
        if any(value is False for value in values):
            return False
    if term.is_Pow:
        base, exponent = term.args
        if exponent == 0:
            return True
        base_nz = _closed_nonzero(base)
        if base_nz is True and exponent.is_finite is not False:
            return True
        if base_nz is False and exponent.is_positive is True:
            return False
        if exponent.is_negative is True and quick_defined_nonzero(base) is True:
            return True
    if term.is_Add and features.has_exp and not exp_indep_inapplicable(term):
        exp_nz = exponential_independence_nonzero(term)
        if exp_nz is True:
            return True
    if (not features.has_float and not features.has_nonfinite
            and features.nodes <= cfg.QUICK_KIND_MAX_NODES):
        kind = number_kind(term)
        if kind is NumberKind.TRANSCENDENTAL:
            return True
    return None


def quick_nonzero(
    term: sp.Expr,
    assumptions=True,
    memo: MutableMapping | None = None,
    *,
    context: OracleMemo | None = None,
) -> bool | None:
    """Return a proof that ``term != 0`` when one is structurally cheap.

    ``True`` means nonzero is proved, ``False`` means zero is proved, and
    ``None`` means this deliberately small engine is inconclusive.
    """
    term = sp.sympify(term)
    assumptions = normalize_assumptions(assumptions)
    key = (term, assumptions)
    cache = context.nonzero_cache if context is not None else memo
    if context is not None and cache is None:
        context.nonzero_cache = {}
        cache = context.nonzero_cache
    if cache is not None and key in cache:
        return cache[key]
    if not term.free_symbols and assumptions is sp.true:
        cached = _closed_nonzero(term)
        if cache is not None:
            cache[key] = cached
        return cached
    features = expression_features(term)
    if features.has_float or features.has_nonfinite:
        return None
    if features.nodes <= 12:
        zero = term.is_zero
        if zero is not None:
            return not bool(zero)
    asked = is_nonzero(term, assumptions)
    if asked is not None:
        return asked
    if term.func is sp.exp and len(term.args) == 1 and term.args[0].is_finite is True:
        return True
    gamma = _gamma_nonzero(term)
    if gamma is not None:
        return gamma
    fact_nonzero = function_nonzero(term, assumptions) if term.func in PROPERTY_FUNCS else None
    if fact_nonzero is True and quick_defined(term, assumptions) is True:
        if cache is not None:
            cache[key] = True
        return True
    if term.is_Add:
        # A sum of nonnegative terms with one strictly positive term is nonzero.
        nonneg = [ask_property(sp.Q.nonnegative, arg, assumptions) for arg in term.args]
        if (
            nonneg
            and all(value is True for value in nonneg)
            and any(is_positive(arg, assumptions) is True for arg in term.args)
        ):
            if cache is not None:
                cache[key] = True
            return True
    if term.is_Mul:
        values = tuple(quick_nonzero(arg, assumptions, memo, context=context) for arg in term.args)
        if all(value is True for value in values):
            return True
        if any(value is False for value in values):
            return False
    if term.is_Pow:
        base, exponent = term.args
        if exponent == 0:
            result = True
            if cache is not None:
                cache[key] = result
            return result
        base_nz = quick_nonzero(base, assumptions, memo, context=context)
        if base_nz is True and exponent.is_finite is not False:
            result = True
            if cache is not None:
                cache[key] = result
            return result
        if exponent.is_negative is True and quick_defined_nonzero(base, assumptions, context=context) is True:
            if cache is not None:
                cache[key] = True
            return True
    if (not term.free_symbols and term.is_Add and features.has_exp
            and not exp_indep_inapplicable(term)):
        exp_nz = exponential_independence_nonzero(term)
        if exp_nz is True:
            if cache is not None:
                cache[key] = True
            return True
    if cache is not None:
        cache[key] = None
    return None


def clear_nonzero_cache() -> None:
    _closed_nonzero.cache_clear()
