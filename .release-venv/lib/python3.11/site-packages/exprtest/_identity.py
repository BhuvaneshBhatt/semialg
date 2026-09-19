"""Small bounded identity normalizers for the latency-sensitive zero oracle."""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from . import _config as cfg
from ._assumptions import is_positive, normalize_assumptions
from ._features import expression_features


def _small(term: sp.Expr, *, ops: int, nodes: int) -> bool:
    feats = expression_features(term)
    return feats.ops <= ops and feats.nodes <= nodes


def _rational_candidate(term: sp.Expr) -> bool:
    """Return whether ``cancel`` is useful and safely bounded for ``term``."""
    if not _small(term, ops=cfg.IDENTITY_RAT_MAX_OPS, nodes=cfg.IDENTITY_RAT_MAX_NODES):
        return False
    if len(term.free_symbols) > cfg.IDENTITY_RAT_MAX_VARS:
        return False
    for node in sp.preorder_traversal(term):
        if node.is_Add and len(node.args) > cfg.IDENTITY_RAT_MAX_ADD_TERMS:
            return False
        if (
            node.is_Pow
            and node.exp.is_Integer
            and abs(int(node.exp)) > cfg.IDENTITY_RAT_MAX_POW_EXP
        ):
            return False
    return True


@lru_cache(maxsize=cfg.NORMAL_FORM_CACHE_SIZE)
def rational_identity_normal_form(term: sp.Expr) -> sp.Expr:
    """Canonicalize a small polynomial/rational expression with ``cancel``.

    Admission is deliberately strict.  ``cancel`` is very effective for the
    common ``P/Q`` identity case but can become expensive on large expanded
    polynomials, so this stage never sees expressions beyond the small limits
    in :mod:`exprtest._config`.
    """
    term = sp.sympify(term)
    if term.has(sp.Float) and not term.free_symbols:
        return term
    if not _rational_candidate(term):
        return term
    try:
        result = sp.cancel(term)
    except (sp.PolynomialError, ValueError, TypeError):
        return term
    if not _small(
        result,
        ops=cfg.IDENTITY_RAT_MAX_RESULT_OPS,
        nodes=cfg.IDENTITY_RAT_MAX_RESULT_NODES,
    ):
        return term
    return result


def _positive_rational_power(base: sp.Expr, exp: sp.Expr) -> sp.Expr:
    """Rewrite powers of positive rationals into prime-power products."""
    if not base.is_Rational or base <= 0 or base == 1:
        return sp.Pow(base, exp)
    num, den = int(base.p), int(base.q)
    factors: dict[int, int] = {}
    for prime, power in sp.factorint(num).items():
        factors[int(prime)] = factors.get(int(prime), 0) + int(power)
    for prime, power in sp.factorint(den).items():
        factors[int(prime)] = factors.get(int(prime), 0) - int(power)
    parts = [
        sp.Pow(sp.Integer(prime), sp.Mul(power, exp), evaluate=False)
        for prime, power in sorted(factors.items())
        if power
    ]
    return sp.Mul(*parts) if parts else sp.Integer(1)


def _identity_node(term: sp.Expr, assumptions=sp.true) -> sp.Expr:
    """Apply a whitelist of branch-safe local power/trig/log identities."""
    if not term.args:
        return term
    args = tuple(_identity_node(arg, assumptions) for arg in term.args)
    try:
        cur = term.func(*args)
    except (TypeError, ValueError):
        cur = term

    # Principal I**z uses log(I) = I*pi/2 exactly.
    if cur.is_Pow and cur.base == sp.I:
        return sp.exp(sp.I * sp.pi * cur.exp / 2)

    # Positive rational bases have an unambiguous real logarithm, so exponent
    # laws are valid for arbitrary complex exponents.
    if cur.is_Pow and cur.base.is_Rational and cur.base > 0:
        return _positive_rational_power(cur.base, cur.exp)

    # (a**b)**c = a**(b*c) is branch-safe for positive a.
    if (
        cur.is_Pow
        and cur.base.is_Pow
        and is_positive(cur.base.base, assumptions) is True
    ):
        return sp.Pow(cur.base.base, cur.base.exp * cur.exp)

    # log(a**x) = x*log(a) is safe for positive real a when x is real.
    if cur.func is sp.log and len(cur.args) == 1:
        arg = cur.args[0]
        if (
            arg.is_Pow
            and is_positive(arg.base, assumptions) is True
            and arg.exp.is_real is True
        ):
            return arg.exp * sp.log(arg.base)

    if cur.func is sp.sin and len(cur.args) == 1:
        arg = cur.args[0]
        # Exact phase shift, restricted to syntactic +/- pi/2.
        rest = sp.expand(arg - sp.pi / 2)
        if not rest.has(sp.pi):
            return sp.cos(rest)
        # Double-angle identity in the direction useful for cancellation.
        coeff, inner = arg.as_coeff_Mul(rational=True)
        if coeff == 2:
            return 2 * sp.sin(inner) * sp.cos(inner)

    if cur.func is sp.cos and len(cur.args) == 1:
        arg = cur.args[0]
        rest = sp.expand(arg - sp.pi / 2)
        if not rest.has(sp.pi):
            return -sp.sin(rest)
        coeff, inner = arg.as_coeff_Mul(rational=True)
        if coeff == 2:
            return sp.cos(inner) ** 2 - sp.sin(inner) ** 2

    if cur.func is sp.tan and len(cur.args) == 1:
        arg = cur.args[0]
        # tan is replaced only when cos(arg) is structurally known nonzero.
        if sp.cos(arg).is_zero is False:
            return sp.sin(arg) / sp.cos(arg)

    return cur


@lru_cache(maxsize=cfg.NORMAL_FORM_CACHE_SIZE)
def elementary_identity_normal_form(term: sp.Expr, assumptions=sp.true) -> sp.Expr:
    """Apply small exact identities, then bounded rational canonicalization."""
    term = sp.sympify(term)
    assumptions = normalize_assumptions(assumptions)
    if term.has(sp.Float) and not term.free_symbols:
        return term
    if not _small(
        term, ops=cfg.IDENTITY_ELEM_MAX_OPS, nodes=cfg.IDENTITY_ELEM_MAX_NODES
    ):
        return term
    rewritten = _identity_node(term, assumptions)

    # Canonicalize both the whole result and small algebraic arguments of
    # elementary functions.  This catches exp((1+r)^-1+r/(1+r)) cheaply.
    def canon(node: sp.Expr) -> sp.Expr:
        if not node.args:
            return node
        vals = tuple(canon(arg) for arg in node.args)
        try:
            rebuilt = node.func(*vals)
        except (TypeError, ValueError):
            rebuilt = node
        return rational_identity_normal_form(rebuilt)

    return canon(rewritten)


def clear_identity_cache() -> None:
    rational_identity_normal_form.cache_clear()
    elementary_identity_normal_form.cache_clear()
