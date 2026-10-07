"""Global metric regularity: Lipschitz, uniform, and Hölder continuity."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import Boundedness, Continuity, DetailedPropertyResult, TruthValue
from ..singularities import continuous
from .boundedness import boundedness
from .common import normalize_domain_set


@dataclass(frozen=True)
class LipschitzResult:
    value: TruthValue
    constant: sp.Expr | None = None
    optimal: TruthValue = TruthValue.UNKNOWN
    domain: object = S.Reals
    source: str | None = None
    evidence: tuple[object, ...] = ()


@dataclass(frozen=True)
class HolderRegularityResult:
    exponent: sp.Expr | None
    constant: sp.Expr | None = None
    optimal: TruthValue = TruthValue.UNKNOWN
    domain: object = S.Reals
    source: str | None = None
    evidence: tuple[object, ...] = ()


def _domain_set(domain, var):
    try:
        return normalize_domain_set(domain, var)
    except (TypeError, ValueError, NotImplementedError):
        return None


def _compact(domain, var):
    d = _domain_set(domain, var)
    return bool(
        d is not None
        and d.is_subset(S.Reals) is True
        and getattr(d, "is_closed", None) is True
        and d.inf.is_finite is True
        and d.sup.is_finite is True
    )


def _derivative_bound(expr, var, domain):
    derivative = sp.diff(expr, var)
    if not derivative.has(var):
        return sp.Abs(derivative), True
    if derivative.func in (sp.sin, sp.cos) or (
        derivative.is_Mul and any(a.func in (sp.sin, sp.cos) for a in derivative.args)
    ):
        coefficient = derivative
        for atom in derivative.atoms(sp.sin, sp.cos):
            coefficient = coefficient / atom
        if not coefficient.has(var):
            return sp.Abs(sp.simplify(coefficient)), True
    d = _domain_set(domain, var)
    if d is None:
        return None, False
    try:
        upper = sp.calculus.util.function_range(sp.Abs(derivative), var, d).sup
        if upper != S.Infinity:
            return sp.simplify(upper), True
    except SYMBOLIC_QUERY_ERRORS:
        pass
    try:
        if boundedness(sp.Abs(derivative), var, domain=domain) is Boundedness.BOUNDED:
            return None, True
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return None, False


def lipschitz(expr, var, *, domain=S.Reals, return_result=False):
    """Certify Lipschitz regularity and return a certified constant when known."""
    expr, var = sp.sympify(expr), sp.sympify(var)
    if not expr.has(var):
        result = LipschitzResult(
            TruthValue.TRUE, S.Zero, TruthValue.TRUE, domain, "constant"
        )
    elif continuous(expr, var, domain=domain) is Continuity.DISCONTINUOUS:
        result = LipschitzResult(
            TruthValue.FALSE, domain=domain, source="discontinuity"
        )
    else:
        constant, ok = _derivative_bound(expr, var, domain)
        result = (
            LipschitzResult(
                TruthValue.TRUE,
                constant,
                TruthValue.UNKNOWN,
                domain,
                "bounded-derivative",
            )
            if ok
            else LipschitzResult(TruthValue.UNKNOWN, domain=domain, source="unresolved")
        )
    if return_result:
        return result
    return (
        result.constant
        if result.value is TruthValue.TRUE and result.constant is not None
        else result.value
    )


def uniformly_continuous(expr, var, *, domain=S.Reals, return_result=False):
    """Certify uniform continuity on the requested domain."""
    expr, var = sp.sympify(expr), sp.sympify(var)
    if domain == S.Reals and expr.is_polynomial(var) and sp.degree(expr, var) > 1:
        continuity = Continuity.CONTINUOUS
        value, source, evidence = TruthValue.FALSE, "nonlinear-polynomial-on-reals", ()
    else:
        continuity = continuous(expr, var, domain=domain)
        value = source = evidence = None
    if value is not None:
        pass
    elif continuity is Continuity.DISCONTINUOUS:
        value, source, evidence = TruthValue.FALSE, "discontinuity", ()
    elif _compact(domain, var) and continuity is Continuity.CONTINUOUS:
        value, source, evidence = TruthValue.TRUE, "heine-cantor", ()
    else:
        lip = lipschitz(expr, var, domain=domain, return_result=True)
        if lip.value is TruthValue.TRUE:
            value, source, evidence = TruthValue.TRUE, "lipschitz", (lip,)
        else:
            value, source, evidence = TruthValue.UNKNOWN, "unresolved", ()
    return (
        DetailedPropertyResult(value, domain=domain, source=source, evidence=evidence)
        if return_result
        else value
    )


def _power_holder(expr, var, domain):
    if isinstance(expr, sp.Pow):
        base, exponent = expr.args
    elif expr.func is sp.Abs and isinstance(expr.args[0], sp.Pow):
        base, exponent = expr.args[0].args
    else:
        return None
    if base not in (var, sp.Abs(var)) or exponent.is_real is not True:
        return None
    if (
        sp.ask(sp.Q.positive(exponent)) is not True
        or sp.ask(sp.Q.lt(exponent, 1)) is not True
    ):
        return None
    d = _domain_set(domain, var)
    if d is None or d.contains(0) not in (True, S.true):
        return None
    return exponent


def holder_regularity(expr, var, *, domain=S.Reals, return_result=False):
    """Return the strongest certified Hölder exponent in ``(0, 1]``."""
    expr, var = sp.sympify(expr), sp.sympify(var)
    exponent = _power_holder(expr, var, domain)
    if exponent is not None:
        result = HolderRegularityResult(
            exponent, S.One, TruthValue.TRUE, domain, "power-cusp"
        )
    else:
        lip = lipschitz(expr, var, domain=domain, return_result=True)
        result = (
            HolderRegularityResult(
                S.One, lip.constant, lip.optimal, domain, "lipschitz", (lip,)
            )
            if lip.value is TruthValue.TRUE
            else HolderRegularityResult(None, domain=domain, source="unresolved")
        )
    return result if return_result else result.exponent
