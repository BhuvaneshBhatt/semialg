"""Executable guarded special-function facts already represented by the catalogue."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .results import Convexity, Monotonicity, Sign


@dataclass(frozen=True)
class GuardedTheorem:
    expression: object
    variable: object
    property: str
    value: object
    guard: object
    theorem_id: str
    domain: object

    def applicable(self, assumptions=sp.S.true):
        # Deliberately conservative: only syntactic truth or SymPy-provable implication.
        g = sp.sympify(self.guard)
        a = sp.sympify(assumptions)
        if g is sp.S.true:
            return True
        try:
            bad = sp.simplify_logic(sp.And(a, sp.Not(g)))
            return bad is sp.S.false
        except (TypeError, ValueError, NotImplementedError):
            return False


def executable_special_theorems(expr, var):
    f = sp.sympify(expr)
    x = sp.sympify(var)
    out = []
    pos = x > 0
    if f == sp.erf(x):
        out += [
            GuardedTheorem(
                f,
                x,
                "monotonicity",
                Monotonicity.STRICTLY_INCREASING,
                sp.S.true,
                "DLMF-7.10-erf-derivative",
                sp.S.Reals,
            ),
            GuardedTheorem(
                f,
                x,
                "range",
                sp.Interval.open(-1, 1),
                sp.S.true,
                "DLMF-7.2-erf-range",
                sp.S.Reals,
            ),
        ]
    if f == sp.erfc(x):
        out += [
            GuardedTheorem(
                f,
                x,
                "monotonicity",
                Monotonicity.STRICTLY_DECREASING,
                sp.S.true,
                "DLMF-7.10-erfc-derivative",
                sp.S.Reals,
            ),
            GuardedTheorem(
                f,
                x,
                "sign",
                Sign.POSITIVE,
                sp.S.true,
                "DLMF-7.2-erfc-positive",
                sp.S.Reals,
            ),
        ]
    if f.func is sp.polygamma and len(f.args) == 2 and f.args[1] == x:
        n = f.args[0]
        integer = sp.Contains(n, sp.S.Integers)
        nonneg = n >= 0
        # (-1)^(n+1) psi^(n)(x)>0 for n>=1; n=0 digamma has a zero.
        if n.is_integer is True and n.is_nonnegative is True and n != 0:
            val = Sign.POSITIVE if int(n) % 2 else Sign.NEGATIVE
            out.append(
                GuardedTheorem(
                    f,
                    x,
                    "sign",
                    val,
                    pos,
                    f"DLMF-5.15-polygamma-{n}",
                    sp.Interval.open(0, sp.oo),
                )
            )
        elif n.is_number is False:
            out.append(
                GuardedTheorem(
                    f,
                    x,
                    "sign",
                    Sign.UNKNOWN,
                    sp.And(pos, integer, nonneg, n >= 1),
                    "DLMF-5.15-polygamma-guarded",
                    sp.Interval.open(0, sp.oo),
                )
            )
    if f.func in (sp.besseli, sp.besselk) and len(f.args) == 2 and f.args[1] == x:
        nu = f.args[0]
        guard = sp.And(pos, sp.Contains(nu, sp.S.Reals), nu >= 0)
        mono = (
            Monotonicity.STRICTLY_INCREASING
            if f.func is sp.besseli
            else Monotonicity.STRICTLY_DECREASING
        )
        out += [
            GuardedTheorem(
                f,
                x,
                "sign",
                Sign.POSITIVE,
                guard,
                "DLMF-10.37-positive",
                sp.Interval.open(0, sp.oo),
            ),
            GuardedTheorem(
                f,
                x,
                "monotonicity",
                mono,
                guard,
                "DLMF-10.37-monotonicity",
                sp.Interval.open(0, sp.oo),
            ),
        ]
    if f == sp.loggamma(x):
        out.append(
            GuardedTheorem(
                f,
                x,
                "convexity",
                Convexity.STRICTLY_CONVEX,
                pos,
                "DLMF-5-loggamma-convex",
                sp.Interval.open(0, sp.oo),
            )
        )
    return tuple(out)
