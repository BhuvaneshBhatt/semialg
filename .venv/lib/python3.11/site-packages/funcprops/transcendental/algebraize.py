from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..real_analysis.common import is_semialgebraic_expression


@dataclass(frozen=True)
class AlgebraizationResult:
    expression: sp.Expr
    substitutions: tuple[tuple[sp.Symbol, sp.Expr], ...]
    constraints: sp.Basic = S.true
    exact: bool = True
    method: str = "none"
    source_domain: object = S.Reals
    relations: tuple[sp.Basic, ...] = ()
    inverse_substitutions: tuple[tuple[sp.Symbol, sp.Expr], ...] = ()

    @property
    def changed(self):
        return bool(self.substitutions)


def _linear_coeff(expr, var):
    try:
        p = sp.Poly(sp.expand(expr), var)
        if p.degree() <= 1:
            return p.nth(1), p.nth(0)
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return None


def _exp_generator(expr, var):
    exps = list(expr.atoms(sp.exp))
    if not exps:
        return None
    coeffs = []
    for e in exps:
        lc = _linear_coeff(e.args[0], var)
        if lc is None:
            return None
        a, b = lc
        if a == 0 or a.is_rational is not True or b.is_real is not True:
            return None
        coeffs.append((e, a, b))
    slopes = [sp.Rational(a) for _, a, _ in coeffs]
    from math import gcd, lcm

    L = 1
    for a in slopes:
        L = lcm(L, int(a.q))
    ints = [int(a * L) for a in slopes]
    g = abs(ints[0])
    for n in ints[1:]:
        g = gcd(g, abs(n))
    slope = sp.Rational(g, L)
    return coeffs, slope if slope else None


def algebraize_exponentials(expr, var):
    expr = sp.sympify(expr)
    data = _exp_generator(expr, var)
    if data is None:
        return AlgebraizationResult(expr, ())
    coeffs, slope = data
    u = sp.Dummy("u", real=True)
    repl = {}
    for e, a, b in coeffs:
        power = sp.simplify(a / slope)
        if power.is_Integer is not True:
            return AlgebraizationResult(expr, ())
        repl[e] = sp.exp(b) * u**power
    new = sp.cancel(expr.xreplace(repl))
    if not is_semialgebraic_expression(new, u):
        return AlgebraizationResult(expr, ())
    return AlgebraizationResult(
        new,
        ((u, sp.exp(slope * var)),),
        u > 0,
        True,
        "exponential_substitution",
        S.Reals,
        (),
        ((var, sp.log(u) / slope),),
    )


def rewrite_hyperbolic_to_exp(expr):
    expr = sp.sympify(expr)
    if expr.has(sp.sinh, sp.cosh, sp.tanh, sp.coth, sp.sech, sp.csch):
        try:
            return sp.cancel(expr.rewrite(sp.exp))
        except SYMBOLIC_QUERY_ERRORS:
            return expr
    return expr


def algebraize_trigonometric(expr, var, *, source_domain=None):
    """Tangent-half-angle algebraization on a single branch.

    Exact on ``(-pi, pi)`` (or a supplied subset thereof):
    ``u=tan(x/2)``, ``sin(x)=2u/(1+u^2)``, ``cos(x)=(1-u^2)/(1+u^2)``.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if not expr.has(sp.sin(var), sp.cos(var)):
        return AlgebraizationResult(expr, ())
    # Only direct sin(x), cos(x); affine phases are handled by periodic reduction first.
    bad = [a for a in expr.atoms(sp.sin, sp.cos) if a.args[0] != var]
    if bad:
        return AlgebraizationResult(expr, ())
    u = sp.Dummy("t", real=True)
    new = expr.xreplace(
        {sp.sin(var): 2 * u / (1 + u**2), sp.cos(var): (1 - u**2) / (1 + u**2)}
    )
    new = sp.cancel(new)
    # rational in u is semialgebraic away from its denominator; 1+u^2 is positive.
    dom = (
        sp.Interval.open(-sp.pi, sp.pi)
        if source_domain is None
        else sp.Intersection(sp.sympify(source_domain), sp.Interval.open(-sp.pi, sp.pi))
    )
    return AlgebraizationResult(
        new,
        ((u, sp.tan(var / 2)),),
        S.true,
        True,
        "tangent_half_angle",
        dom,
        (),
        ((var, 2 * sp.atan(u)),),
    )


def algebraize_trig_circle(expr, var):
    """Algebraize direct sin/cos dependence over a full phase using s^2+c^2=1.

    This transform is exact for existential range/sign questions over a full
    period, but callers must respect the returned relation.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if not (expr.has(sp.sin(var)) or expr.has(sp.cos(var))):
        return AlgebraizationResult(expr, ())
    if any(a.args[0] != var for a in expr.atoms(sp.sin, sp.cos)):
        return AlgebraizationResult(expr, ())
    s = sp.Dummy("s", real=True)
    c = sp.Dummy("c", real=True)
    new = expr.xreplace({sp.sin(var): s, sp.cos(var): c})
    return AlgebraizationResult(
        new,
        ((s, sp.sin(var)), (c, sp.cos(var))),
        sp.Eq(s**2 + c**2, 1),
        True,
        "trig_circle",
        sp.Interval.Ropen(0, 2 * sp.pi),
        (sp.Eq(s**2 + c**2, 1),),
    )


def algebraize_common_subexprs(expr, var):
    expr = sp.sympify(expr)
    counts = {}
    for node in sp.preorder_traversal(expr):
        if isinstance(node, sp.Function) and node.has(var):
            counts[node] = counts.get(node, 0) + 1
    repeated = [e for e, n in counts.items() if n > 1]
    if not repeated:
        return AlgebraizationResult(expr, ())
    repl = {
        e: sp.Dummy(f"t{i}", real=True)
        for i, e in enumerate(sorted(repeated, key=sp.default_sort_key), 1)
    }
    return AlgebraizationResult(
        expr.xreplace(repl),
        tuple((v, k) for k, v in repl.items()),
        S.true,
        False,
        "common_subexpression",
    )


def algebraize(expr, var, *, domain=None, allow_trig=False):
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if is_semialgebraic_expression(expr, var):
        return AlgebraizationResult(expr, ())
    rewritten = rewrite_hyperbolic_to_exp(expr)
    r = algebraize_exponentials(rewritten, var)
    if r.changed:
        return r
    if allow_trig:
        r = algebraize_trigonometric(expr, var, source_domain=domain)
        if r.changed:
            return r
    return AlgebraizationResult(expr, ())
