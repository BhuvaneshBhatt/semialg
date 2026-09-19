"""Cached one-pass structural fingerprints for oracle routing."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import sympy as sp

from . import _config as cfg


@dataclass(frozen=True)
class ExprFeatures:
    """Cheap structural facts used to route fast zero-test stages."""

    nodes: int
    depth: int
    ops: int
    generators: int
    roots: int
    logs: int
    max_pow_bits: int
    max_int_bits: int
    has_log: bool
    has_exp: bool
    has_trig: bool
    has_pi: bool
    has_special: bool
    has_radical: bool
    has_algnum: bool
    has_rootof: bool
    has_float: bool
    has_nonfinite: bool
    has_rewrite: bool = False
    has_domain_hazard: bool = False

    @property
    def has_log_exp(self) -> bool:
        return self.has_log or self.has_exp

    @property
    def has_cyclotomic_shape(self) -> bool:
        return self.has_pi or self.has_exp or self.has_trig

    @property
    def has_tower_shape(self) -> bool:
        return self.has_radical or self.has_algnum or self.has_rootof


_SPECIAL_FUNCS = {
    sp.gamma,
    sp.zeta,
    sp.factorial,
    sp.binomial,
    sp.beta,
    sp.factorial2,
    sp.rf,
    sp.ff,
    sp.harmonic,
}
_TRIG_FUNCS = {sp.sin, sp.cos, sp.tan}
_REWRITE_FUNCS = {
    sp.Abs,
    sp.sign,
    sp.conjugate,
    sp.re,
    sp.im,
    sp.sinh,
    sp.cosh,
    sp.tanh,
}


@lru_cache(maxsize=cfg.EXACT_FEATURE_CACHE_SIZE)
def expression_features(term: sp.Expr) -> ExprFeatures:
    """Return a cached structural fingerprint after one explicit DAG walk."""
    term = sp.sympify(term)
    nodes = 0
    depth = 0
    ops = 0
    generators = 0
    roots = 0
    logs = 0
    max_pow_bits = 0
    max_int_bits = 0
    has_log = has_exp = has_trig = has_pi = False
    has_special = has_radical = has_algnum = has_rootof = False
    has_float = has_nonfinite = has_rewrite = has_domain_hazard = False
    stack = [(term, 1)]
    while stack:
        cur, dep = stack.pop()
        nodes += 1
        depth = max(depth, dep)
        if cur.is_Add or cur.is_Mul:
            ops += max(1, len(cur.args) - 1)
        elif cur.args:
            ops += 1
        if cur.is_Integer:
            max_int_bits = max(max_int_bits, abs(int(cur)).bit_length())
        if cur.is_Pow:
            if cur.exp.is_Integer:
                max_pow_bits = max(max_pow_bits, abs(int(cur.exp)).bit_length())
            if cur.exp.is_integer is not True or cur.exp.is_negative is True:
                has_domain_hazard = True
        if cur.func is sp.log:
            has_log = True
            has_domain_hazard = True
            logs += 1
        elif cur.func is sp.exp:
            has_exp = True
        elif cur.func in _TRIG_FUNCS:
            has_trig = True
            if cur.func is sp.tan:
                has_domain_hazard = True
        if cur is sp.pi:
            has_pi = True
        if cur.func in _SPECIAL_FUNCS:
            has_special = True
            if cur.func in {sp.gamma, sp.factorial}:
                has_domain_hazard = True
        if cur.func in _REWRITE_FUNCS:
            has_rewrite = True
            if cur.func is sp.tanh:
                has_domain_hazard = True
        elif cur.func in {sp.atan, sp.loggamma}:
            has_domain_hazard = True
        is_root = bool(getattr(cur, "is_CRootOf", False))
        is_radical = bool(cur.is_Pow and cur.exp.is_Rational and not cur.exp.is_Integer)
        is_algnum = isinstance(cur, sp.AlgebraicNumber)
        if is_root or is_radical or is_algnum:
            generators += 1
        if is_root:
            roots += 1
            has_rootof = True
        if is_radical:
            has_radical = True
        if is_algnum:
            has_algnum = True
        if cur.is_Float:
            has_float = True
        if cur in (sp.nan, sp.zoo, sp.oo, -sp.oo):
            has_nonfinite = True
            has_domain_hazard = True
        stack.extend((arg, dep + 1) for arg in cur.args)
    return ExprFeatures(
        nodes,
        depth,
        ops,
        generators,
        roots,
        logs,
        max_pow_bits,
        max_int_bits,
        has_log,
        has_exp,
        has_trig,
        has_pi,
        has_special,
        has_radical,
        has_algnum,
        has_rootof,
        has_float,
        has_nonfinite,
        has_rewrite,
        has_domain_hazard,
    )


def clear_feature_cache() -> None:
    expression_features.cache_clear()
