"""Verified affine and nonlinear interval Krawczyk certificates."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


def _ival(expr, xs, box):
    e = sp.sympify(expr)
    if not e.has(*xs):
        return (e, e)
    if e in xs:
        interval = box[xs.index(e)]
        return (sp.sympify(interval.start), sp.sympify(interval.end))
    if e.is_Add:
        parts = [_ival(a, xs, box) for a in e.args]
        return (
            sp.simplify(sum(a for a, b in parts)),
            sp.simplify(sum(b for a, b in parts)),
        )
    if e.is_Mul:
        lo = hi = sp.Integer(1)
        for a in e.args:
            al, ah = _ival(a, xs, box)
            vals = [lo * al, lo * ah, hi * al, hi * ah]
            try:
                lo, hi = min(vals), max(vals)
            except TypeError:
                return None
        return (sp.simplify(lo), sp.simplify(hi))
    if e.is_Pow and e.exp.is_Integer and e.exp >= 0:
        base = _ival(e.base, xs, box)
        if base is None:
            return None
        lower, upper = base
        n = int(e.exp)
        vals = [lower**n, upper**n]
        if n % 2 == 0 and lower <= 0 <= upper:
            vals.append(sp.Integer(0))
        try:
            return (sp.simplify(min(vals)), sp.simplify(max(vals)))
        except TypeError:
            return None
    return None


def _add(A, B):
    return (sp.simplify(A[0] + B[0]), sp.simplify(A[1] + B[1]))


def _mul(A, B):
    vals = [A[0] * B[0], A[0] * B[1], A[1] * B[0], A[1] * B[1]]
    try:
        return (sp.simplify(min(vals)), sp.simplify(max(vals)))
    except TypeError:
        return None


def _krawczyk(expressions, xs, box):
    F = sp.Matrix(expressions)
    J = F.jacobian(xs)
    x0 = sp.Matrix([(interval.start + interval.end) / 2 for interval in box])
    sub = dict(zip(xs, x0, strict=False))
    J0 = J.subs(sub)
    try:
        C = J0.inv()
    except (ValueError, ZeroDivisionError, NotImplementedError):
        return None
    f0 = F.subs(sub)
    base = x0 - C * f0
    n = len(xs)
    JI = [[_ival(J[i, j], xs, box) for j in range(n)] for i in range(n)]
    if any(v is None for row in JI for v in row):
        return None
    M = []
    for i in range(n):
        row = []
        for j in range(n):
            s = (sp.Integer(0), sp.Integer(0))
            for k in range(n):
                q = _mul((C[i, k], C[i, k]), JI[k][j])
                if q is None:
                    return None
                s = _add(s, q)
            delta = (sp.Integer(1) if i == j else sp.Integer(0),) * 2
            row.append((sp.simplify(delta[0] - s[1]), sp.simplify(delta[1] - s[0])))
        M.append(row)
    dx = [
        (interval.start - x0[j], interval.end - x0[j]) for j, interval in enumerate(box)
    ]
    K = []
    for i in range(n):
        s = (base[i], base[i])
        for j in range(n):
            q = _mul(M[i][j], dx[j])
            if q is None:
                return None
            s = _add(s, q)
        K.append(s)
    return tuple(K), tuple(x0)


@dataclass(frozen=True)
class KrawczykBoxCertificate:
    expressions: tuple
    variables: tuple
    box: tuple
    root: tuple
    unique: bool
    source: str = "interval-krawczyk"
    enclosure: tuple = ()

    def verify(self):
        if self.source == "affine-krawczyk":
            F = sp.Matrix(self.expressions)
            J = F.jacobian(self.variables)
            if any(e.free_symbols.intersection(self.variables) for e in J):
                return False
            return self.unique and all(
                interval.contains(v) is sp.S.true
                for interval, v in zip(self.box, self.root, strict=False)
            )
        if self.source != "interval-krawczyk" or not self.unique:
            return False
        out = _krawczyk(self.expressions, self.variables, self.box)
        if out is None:
            return False
        K, _x0 = out
        if K != self.enclosure:
            return False
        return all(
            (interval.start < k[0]) is sp.S.true and (k[1] < interval.end) is sp.S.true
            for interval, k in zip(self.box, K, strict=False)
        )


def krawczyk_box(expressions, variables, box):
    xs = tuple(map(sp.sympify, variables))
    box = tuple(box)
    F = tuple(map(sp.sympify, expressions))
    if (
        len(F) != len(xs)
        or len(box) != len(xs)
        or not all(
            isinstance(interval, sp.Interval)
            and not interval.left_open
            and not interval.right_open
            for interval in box
        )
    ):
        return None
    out = _krawczyk(F, xs, box)
    if out is None:
        return None
    K, x0 = out
    if not all(
        (interval.start < k[0]) is sp.S.true and (k[1] < interval.end) is sp.S.true
        for interval, k in zip(box, K, strict=False)
    ):
        return None
    c = KrawczykBoxCertificate(F, xs, box, x0, True, "interval-krawczyk", K)
    return c if c.verify() else None
