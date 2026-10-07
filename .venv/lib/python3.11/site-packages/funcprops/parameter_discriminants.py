"""Certified transition loci for parameterized univariate families."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .exactness import certify_identity_equal


@dataclass(frozen=True)
class TranscendentalDiscriminant:
    family: str
    condition: sp.Expr
    transition: sp.Expr
    source: str


@dataclass(frozen=True)
class ParameterTransitionLocus:
    expression: sp.Expr
    variable: sp.Symbol
    parameters: tuple
    locus: sp.Expr
    complete: bool
    source: str

    def verify(self):
        if self.source != "polynomial-resultant":
            return False
        try:
            p = sp.Poly(self.expression, self.variable)
            r = sp.factor(
                sp.resultant(
                    p.as_expr(), sp.diff(p.as_expr(), self.variable), self.variable
                )
            )
        except sp.PolynomialError:
            return False
        return certify_identity_equal(r, self.locus) is True and self.complete


def polynomial_transition_locus(expr, var, parameters=None):
    f = sp.sympify(expr)
    x = sp.sympify(var)
    try:
        p = sp.Poly(f, x)
    except sp.PolynomialError:
        return None
    pars = tuple(
        sorted(
            (f.free_symbols - {x})
            if parameters is None
            else map(sp.sympify, parameters),
            key=str,
        )
    )
    if not pars:
        return None
    r = sp.factor(sp.resultant(p.as_expr(), sp.diff(p.as_expr(), x), x))
    c = ParameterTransitionLocus(f, x, pars, r, True, "polynomial-resultant")
    return c if c.verify() else None


def exp_linear_discriminant(expr, var, a, b):
    x = sp.sympify(var)
    a = sp.sympify(a)
    b = sp.sympify(b)
    f = sp.sympify(expr)
    if certify_identity_equal(f, sp.exp(x) + a * x + b) is not True:
        return None
    D = sp.simplify(b + a * (sp.log(-a) - 1))
    return TranscendentalDiscriminant(
        "exp-linear", a < 0, D, "critical-point-elimination"
    )


def exp_linear_root_strata(expr, var, a, b):
    d = exp_linear_discriminant(expr, var, a, b)
    if d is None:
        return None
    D = d.transition
    return (
        (2, sp.And(a < 0, D < 0)),
        (1, sp.And(a < 0, sp.Eq(D, 0))),
        (0, sp.And(a < 0, D > 0)),
        (1, a > 0),
        (1, sp.And(sp.Eq(a, 0), b < 0)),
        (0, sp.And(sp.Eq(a, 0), b >= 0)),
    )
