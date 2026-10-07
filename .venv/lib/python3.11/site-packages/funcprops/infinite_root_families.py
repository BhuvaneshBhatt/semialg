"""Structured complete infinite real-root families and ZeroGeometry adapters."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .exactness import certify_identity_equal


@dataclass(frozen=True)
class IndexedRootFamily:
    expression: sp.Expr
    variable: sp.Symbol
    formula: sp.Expr
    index: sp.Symbol
    index_domain: object
    simple: bool
    complete: bool
    source: str
    multiplicity: int = 1

    def verify(self):
        if not self.complete or self.multiplicity < 1:
            return False
        x = self.variable
        if certify_identity_equal(self.expression.subs(x, self.formula), 0) is not True:
            return False
        d = sp.diff(self.expression, x).subs(x, self.formula)
        if self.simple and sp.simplify(d).is_zero is not False:
            return False
        return True

    def to_zero_geometry(self):
        if not self.verify():
            return None
        from .transcendental.zero_geometry import ZeroGeometryCertificate

        return ZeroGeometryCertificate(
            self.source,
            sp.S.Reals,
            sp.oo,
            sp.S.true,
            self.simple,
            (),
            (),
            (),
            0,
            self.complete,
            "structured-indexed-family",
            notes=(self.formula, self.index, self.index_domain, self.multiplicity),
            real_exhaustive=True,
        )


def _affine_arg(f, x):
    arg = sp.expand(f.args[0])
    a = sp.simplify(sp.diff(arg, x))
    b = sp.simplify(arg - a * x)
    if a.has(x) or a.is_zero is not False or a.is_real is not True:
        return None
    return a, b


def affine_trig_root_family(expr, var):
    f = sp.sympify(expr)
    x = sp.sympify(var)
    n = sp.Symbol("n", integer=True)
    if f.func not in (sp.sin, sp.cos, sp.tan):
        return None
    ab = _affine_arg(f, x)
    if ab is None:
        return None
    a, b = ab
    phase = 0 if f.func in (sp.sin, sp.tan) else sp.pi / 2
    formula = sp.simplify((sp.pi * n + phase - b) / a)
    source = {
        sp.sin: "affine-sine-zeros",
        sp.cos: "affine-cosine-zeros",
        sp.tan: "affine-tangent-zeros",
    }[f.func]
    fam = IndexedRootFamily(f, x, formula, n, sp.S.Integers, True, True, source)
    return fam if fam.verify() else None


def affine_sine_root_family(expr, var):
    fam = affine_trig_root_family(expr, var)
    return fam if fam and fam.expression.func is sp.sin else None


def structured_zero_geometry(expr, var):
    fam = affine_trig_root_family(expr, var)
    return None if fam is None else fam.to_zero_geometry()
