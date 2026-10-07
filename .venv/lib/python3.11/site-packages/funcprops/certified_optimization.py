"""Assembly of replayable univariate global-minimum certificates."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .exactness import certify_identity_equal, certify_zero


@dataclass(frozen=True)
class GlobalMinimumCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    point: sp.Expr
    value: sp.Expr
    unique: bool
    source: str
    domain: object = sp.S.Reals
    coercive: bool = False

    def verify(self):
        x = self.variable
        f = self.expression
        p = self.point
        D = self.domain
        if D.contains(p) is not sp.S.true:
            return False
        if certify_zero(sp.diff(f, x).subs(x, p)) is not True:
            return False
        if certify_identity_equal(f.subs(x, p), self.value) is not True:
            return False
        d2 = sp.simplify(sp.diff(f, x, 2))
        try:
            lo = sp.calculus.util.minimum(d2, x, D)
        except (ValueError, TypeError, NotImplementedError, sp.PolynomialError):
            return False
        convex = lo is not None and lo is not -sp.oo and lo.is_nonnegative is True
        strict = lo is not None and lo is not -sp.oo and lo.is_positive is True
        if not convex:
            return False
        if self.unique and not strict:
            return False
        if self.coercive and D is sp.S.Reals:
            try:
                if (
                    sp.limit(f, x, sp.oo) is not sp.oo
                    or sp.limit(f, x, -sp.oo) is not sp.oo
                ):
                    return False
            except (ValueError, NotImplementedError):
                return False
        return self.source in {
            "convex-stationary",
            "strict-convex-stationary",
            "coercive-convex-stationary",
        }


def certify_global_minimum(expr, var, domain=sp.S.Reals):
    f = sp.sympify(expr)
    x = sp.sympify(var)
    D = domain if isinstance(domain, sp.Set) else sp.S.Reals
    d1 = sp.diff(f, x)
    d2 = sp.simplify(sp.diff(f, x, 2))
    try:
        lo = sp.calculus.util.minimum(d2, x, D)
    except (ValueError, TypeError, NotImplementedError, sp.PolynomialError):
        return None
    if lo is None or lo is -sp.oo or lo.is_nonnegative is not True:
        return None
    roots = sp.solveset(sp.Eq(d1, 0), x, domain=D)
    if not isinstance(roots, sp.FiniteSet) or len(roots) != 1:
        return None
    p = next(iter(roots))
    strict = lo.is_positive is True
    coercive = False
    if D is sp.S.Reals:
        try:
            coercive = (
                sp.limit(f, x, sp.oo) is sp.oo and sp.limit(f, x, -sp.oo) is sp.oo
            )
        except (ValueError, NotImplementedError):
            coercive = False
    src = (
        "coercive-convex-stationary"
        if coercive
        else "strict-convex-stationary"
        if strict
        else "convex-stationary"
    )
    c = GlobalMinimumCertificate(
        f, x, p, sp.simplify(f.subs(x, p)), strict, src, D, coercive
    )
    return c if c.verify() else None
