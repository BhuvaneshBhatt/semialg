"""Replayable eventual-sign certificates for conservative dominance families."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class TailSignCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    direction: int
    threshold: sp.Expr
    sign: int
    source: str

    def verify(self):
        x = self.variable
        f = sp.sympify(self.expression)
        T = sp.sympify(self.threshold)
        if (
            self.direction not in (-1, 1)
            or self.sign not in (-1, 1)
            or T.is_real is False
        ):
            return False
        y = sp.Symbol("_tail_y", nonnegative=True)
        xx = self.direction * (sp.Abs(T) + y)
        if self.source in {"polynomial-leading-term", "rational-leading-term"}:
            num, den = sp.fraction(sp.cancel(f))
            try:
                _pn, _pd = sp.Poly(num, x), sp.Poly(den, x)
            except sp.PolynomialError:
                return False
            qn = sp.Poly(sp.expand(num.subs(x, xx)), y)
            qd = sp.Poly(sp.expand(den.subs(x, xx)), y)

            def signed(p, s):
                cs = p.all_coeffs()
                return (
                    all((c * s).is_nonnegative is True for c in cs)
                    and (p.TC() * s).is_positive is True
                )

            # denominator must have a fixed nonzero sign on the tail.
            ds = 1 if signed(qd, 1) else -1 if signed(qd, -1) else 0
            return ds != 0 and signed(qn, self.sign * ds)
        if self.source == "exp-polynomial-dominance":
            terms = sp.Mul.make_args(f)
            ex = [t for t in terms if t.func is sp.exp]
            if len(ex) != 1:
                return False
            a = sp.simplify(sp.diff(ex[0].args[0], x))
            if a.has(x) or (a * self.direction).is_positive is not True:
                return False
            rest = sp.simplify(f / ex[0])
            try:
                sp.Poly(rest, x)
            except sp.PolynomialError:
                return False
            q = sp.Poly(sp.expand(rest.subs(x, xx)), y)
            return (
                all((c * self.sign).is_nonnegative is True for c in q.all_coeffs())
                and (q.TC() * self.sign).is_positive is True
            )
        return False


def certify_tail_sign(expr, var, *, direction=1, threshold=None):
    f = sp.sympify(expr)
    x = sp.sympify(var)
    candidates = (
        [0, 1, 2, 4, 8, 16, 32, 64] if threshold is None else [sp.sympify(threshold)]
    )
    num, den = sp.fraction(sp.cancel(f))
    try:
        pn, pd = sp.Poly(num, x), sp.Poly(den, x)
    except sp.PolynomialError:
        pn = pd = None
    if (
        pn is not None
        and pd is not None
        and pn.LC().is_nonzero is True
        and pd.LC().is_nonzero is True
    ):
        s = (
            1
            if (pn.LC() / pd.LC()).is_positive
            else -1
            if (pn.LC() / pd.LC()).is_negative
            else None
        )
        if s is not None and direction < 0 and (pn.degree() - pd.degree()) % 2:
            s = -s
        if s is not None:
            src = (
                "polynomial-leading-term"
                if pd.degree() == 0
                else "rational-leading-term"
            )
            for T in candidates:
                c = TailSignCertificate(f, x, direction, sp.sympify(T), s, src)
                if c.verify():
                    return c
    # exp(a*x+b)*P(x), where exponential is positive and grows in the requested direction.
    for e in sp.Mul.make_args(f):
        if e.func is sp.exp:
            a = sp.simplify(sp.diff(e.args[0], x))
            rest = sp.simplify(f / e)
            try:
                p = sp.Poly(rest, x)
            except sp.PolynomialError:
                continue
            if (
                a.has(x)
                or (a * direction).is_positive is not True
                or p.LC().is_nonzero is not True
            ):
                continue
            s = 1 if p.LC().is_positive else -1
            if direction < 0 and p.degree() % 2:
                s = -s
            for T in candidates:
                c = TailSignCertificate(
                    f, x, direction, sp.sympify(T), s, "exp-polynomial-dominance"
                )
                if c.verify():
                    return c
    return None
