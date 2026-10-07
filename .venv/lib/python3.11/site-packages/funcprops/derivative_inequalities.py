"""Derivative-chain inequality proofs with replayable endpoint induction."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .certified_inequalities import certify_inequality
from .results import TruthValue


@dataclass(frozen=True)
class DerivativeChainCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    domain: sp.Interval
    depth: int
    endpoint: sp.Expr
    base_certificate: object

    def verify(self):
        if not self.base_certificate.verify():
            return False
        x = self.variable
        f = self.expression
        a = self.endpoint
        # If f^(depth)>0 and all lower derivatives are nonnegative at left endpoint,
        # repeated integration proves f>0 on [a,b] (strict for interior).
        return all(
            sp.simplify(sp.diff(f, x, k).subs(x, a)).is_nonnegative is True
            for k in range(self.depth)
        )


def derivative_chain_positive(expr, var, domain, *, max_depth=6):
    if not isinstance(domain, sp.Interval) or domain.start in (-sp.oo, sp.oo):
        return None
    f = sp.sympify(expr)
    x = sp.sympify(var)
    for d in range(max_depth + 1):
        g = sp.diff(f, x, d)
        cert = certify_inequality(g, 0, x, domain, relation=">", max_depth=5)
        if cert.classification is TruthValue.TRUE:
            out = DerivativeChainCertificate(f, x, domain, d, domain.start, cert)
            if out.verify():
                return out
    return None
