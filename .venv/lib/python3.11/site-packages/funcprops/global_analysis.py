"""Compact-core plus certified-tail decomposition for global inequalities."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .asymptotic_tails import TailSignCertificate, certify_tail_sign
from .certified_inequalities import InequalityCertificate, certify_inequality
from .results import TruthValue


@dataclass(frozen=True)
class GlobalInequalityCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    relation: str
    core: sp.Interval
    core_certificate: InequalityCertificate
    left_tail: TailSignCertificate | None
    right_tail: TailSignCertificate | None
    complete: bool

    def verify(self):
        if not self.core_certificate.verify():
            return False
        if self.left_tail is not None and not self.left_tail.verify():
            return False
        if self.right_tail is not None and not self.right_tail.verify():
            return False
        return self.complete == (
            self.left_tail is not None and self.right_tail is not None
        )


def certify_global_positive(expr, var, *, core_radius=8):
    x = sp.sympify(var)
    f = sp.sympify(expr)
    R = sp.sympify(core_radius)
    core = sp.Interval(-R, R)
    cc = certify_inequality(f, 0, x, core, relation=">", max_depth=10)
    lt = certify_tail_sign(f, x, direction=-1, threshold=-R)
    rt = certify_tail_sign(f, x, direction=1, threshold=R)
    complete = (
        cc.classification is TruthValue.TRUE and lt is not None and rt is not None
    )
    return GlobalInequalityCertificate(f, x, ">", core, cc, lt, rt, complete)
