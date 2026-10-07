"""Argument-principle and Rouché certificates on circles."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy.polys.polyerrors import DomainError

from .exactness import certify_identity_equal


@dataclass(frozen=True)
class RoucheCertificate:
    f: sp.Expr
    g: sp.Expr
    var: sp.Symbol
    center: sp.Expr
    radius: sp.Expr
    count: int
    source: str = "polynomial-circle-rouche"

    def verify(self):
        if self.source != "polynomial-circle-rouche":
            return False
        # Current replayable theorem: dominant monomial a*(z-c)^n against polynomial remainder.
        z = self.var
        sp.Symbol("_w", real=True, nonnegative=True)
        shifted = sp.expand(self.f.subs(z, z + self.center))
        try:
            p = sp.Poly(shifted, z)
        except sp.PolynomialError:
            return False
        n = p.degree()
        a = p.LC()
        rem = sp.expand(shifted - a * z**n)
        if n != self.count or not self.radius.is_positive:
            return False
        if certify_identity_equal(self.g, a * (z - self.center) ** n) is not True:
            return False
        # Triangle-bound sufficient condition on |z|=R.
        bound = (
            sp.Add(*[abs(c) * self.radius**k for (k,), c in sp.Poly(rem, z).terms()])
            if rem != 0
            else 0
        )
        dominant = abs(a) * self.radius**n
        q = sp.simplify(dominant - bound)
        return q.is_positive is True


def rouche_polynomial_circle(expr, var, *, center=0, radius=1):
    z = sp.sympify(var)
    f = sp.sympify(expr)
    c = sp.sympify(center)
    R = sp.sympify(radius)
    shifted = sp.expand(f.subs(z, z + c))
    try:
        p = sp.Poly(shifted, z)
    except sp.PolynomialError:
        return None
    n = p.degree()
    a = p.LC()
    rem = sp.expand(shifted - a * z**n)
    if R.is_positive is not True:
        return None
    terms = sp.Poly(rem, z).terms() if rem != 0 else []
    bound = sp.Add(*[abs(coef) * R**k for (k,), coef in terms]) if terms else 0
    if sp.simplify(abs(a) * R**n - bound).is_positive is True:
        return RoucheCertificate(f, a * (z - c) ** n, z, c, R, n)
    return None


@dataclass(frozen=True)
class ArgumentPrincipleCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    lower_left: sp.Expr
    upper_right: sp.Expr
    count: int
    source: str = "polynomial-argument-principle"

    def verify(self):
        if self.source != "polynomial-argument-principle":
            return False
        try:
            p = sp.Poly(self.expression, self.variable)
            return int(p.count_roots(self.lower_left, self.upper_right)) == self.count
        except (
            sp.PolynomialError,
            DomainError,
            TypeError,
            ValueError,
            NotImplementedError,
        ):
            return False


def argument_principle_rectangle(expr, var, lower_left, upper_right):
    """Exact polynomial zero count in a complex rectangle.

    SymPy's exact polynomial root counter is an argument-principle/Cauchy-index
    computation; the returned certificate independently replays the count.
    """
    try:
        p = sp.Poly(sp.sympify(expr), sp.sympify(var))
        n = int(p.count_roots(sp.sympify(lower_left), sp.sympify(upper_right)))
    except (
        sp.PolynomialError,
        DomainError,
        TypeError,
        ValueError,
        NotImplementedError,
    ):
        return None
    return ArgumentPrincipleCertificate(
        sp.sympify(expr),
        sp.sympify(var),
        sp.sympify(lower_left),
        sp.sympify(upper_right),
        n,
    )


@dataclass(frozen=True)
class RectangleRoucheCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    center: sp.Expr
    lower_left: sp.Expr
    upper_right: sp.Expr
    count: int
    dominant: sp.Expr
    source: str = "polynomial-rectangle-rouche"

    def verify(self):
        if self.source != "polynomial-rectangle-rouche":
            return False
        c = rouche_polynomial_rectangle(
            self.expression, self.variable, self.lower_left, self.upper_right
        )
        return (
            c is not None
            and c.count == self.count
            and certify_identity_equal(c.dominant, self.dominant) is True
        )


def rouche_polynomial_rectangle(expr, var, lower_left, upper_right):
    """Rouché count using a centered dominant monomial on a rectangle boundary."""
    z = sp.sympify(var)
    f = sp.sympify(expr)
    ll = sp.sympify(lower_left)
    ur = sp.sympify(upper_right)
    xmin, xmax, ymin, ymax = map(
        sp.simplify, (sp.re(ll), sp.re(ur), sp.im(ll), sp.im(ur))
    )
    if not all(q.is_real is True for q in (xmin, xmax, ymin, ymax)):
        return None
    hx = sp.simplify((xmax - xmin) / 2)
    hy = sp.simplify((ymax - ymin) / 2)
    if hx.is_positive is not True or hy.is_positive is not True:
        return None
    c = sp.simplify((xmin + xmax) / 2 + sp.I * (ymin + ymax) / 2)
    w = sp.Symbol("_rw")
    try:
        p = sp.Poly(sp.expand(f.subs(z, w + c)), w)
    except sp.PolynomialError:
        return None
    n = p.degree()
    a = p.LC()
    rem = sp.expand(p.as_expr() - a * w**n)
    rmin = sp.Min(hx, hy)
    rmax = sp.sqrt(hx**2 + hy**2)
    terms = sp.Poly(rem, w).terms() if rem != 0 else []
    bound = sp.Add(*[abs(co) * rmax**k for (k,), co in terms]) if terms else 0
    if sp.simplify(abs(a) * rmin**n - bound).is_positive is not True:
        return None
    return RectangleRoucheCertificate(f, z, c, ll, ur, n, a * (z - c) ** n)
