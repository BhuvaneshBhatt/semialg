"""Exact meromorphic zero-minus-pole counts for rational functions on rectangles."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy.polys.polyerrors import DomainError


def _boundary_clear(poly, lower_left, upper_right):
    x1, y1 = sp.re(lower_left), sp.im(lower_left)
    x2, y2 = sp.re(upper_right), sp.im(upper_right)
    corners = (x1 + sp.I * y1, x2 + sp.I * y1, x2 + sp.I * y2, x1 + sp.I * y2)
    try:
        return all(
            int(poly.count_roots(a, b)) == 0
            for a, b in zip(corners, corners[1:] + corners[:1], strict=False)
        )
    except (DomainError, ValueError, TypeError, NotImplementedError):
        return False


@dataclass(frozen=True)
class MeromorphicContourCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    lower_left: sp.Expr
    upper_right: sp.Expr
    zeros: int
    poles: int
    winding: int
    source: str = "rational-argument-principle"

    def verify(self):
        if self.source != "rational-argument-principle":
            return False
        num, den = map(sp.cancel, (sp.fraction(sp.cancel(self.expression))))
        try:
            pn = sp.Poly(num, self.variable)
            pd = sp.Poly(den, self.variable)
            if not _boundary_clear(
                pn, self.lower_left, self.upper_right
            ) or not _boundary_clear(pd, self.lower_left, self.upper_right):
                return False
            z = int(pn.count_roots(self.lower_left, self.upper_right))
            p = int(pd.count_roots(self.lower_left, self.upper_right))
        except (
            sp.PolynomialError,
            DomainError,
            ValueError,
            TypeError,
            NotImplementedError,
        ):
            return False
        return (z, p, z - p) == (self.zeros, self.poles, self.winding)


def meromorphic_rectangle_count(expr, var, lower_left, upper_right):
    f = sp.cancel(sp.sympify(expr))
    x = sp.sympify(var)
    num, den = sp.fraction(f)
    try:
        pn = sp.Poly(num, x)
        pd = sp.Poly(den, x)
        if not _boundary_clear(pn, lower_left, upper_right) or not _boundary_clear(
            pd, lower_left, upper_right
        ):
            return None
        z = int(pn.count_roots(lower_left, upper_right))
        p = int(pd.count_roots(lower_left, upper_right))
    except (
        sp.PolynomialError,
        DomainError,
        ValueError,
        TypeError,
        NotImplementedError,
    ):
        return None
    c = MeromorphicContourCertificate(
        f, x, sp.sympify(lower_left), sp.sympify(upper_right), z, p, z - p
    )
    return c if c.verify() else None
