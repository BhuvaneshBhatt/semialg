"""Certified real inequality engine."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .certified_balls import certified_real_ball
from .results import TruthValue
from .taylor_models import taylor_model


@dataclass(frozen=True)
class InequalityCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    domain: sp.Interval
    relation: str
    classification: TruthValue
    source: str
    evidence: object = None

    def verify(self):
        b = certified_real_ball(self.expression, self.variable, self.domain)
        if b is None:
            return False
        if self.relation == ">":
            return (self.classification is TruthValue.TRUE) == (
                b.lower.is_positive is True
            )
        if self.relation == "<":
            return (self.classification is TruthValue.TRUE) == (
                b.upper.is_negative is True
            )
        if self.relation == ">=":
            return (self.classification is TruthValue.TRUE) == (
                b.lower.is_nonnegative is True
            )
        if self.relation == "<=":
            return (self.classification is TruthValue.TRUE) == (
                b.upper.is_nonpositive is True
            )
        return False


def certify_inequality(lhs, rhs, var, domain, *, relation=">", max_depth=8):
    f = sp.simplify(sp.sympify(lhs) - sp.sympify(rhs))
    x = sp.sympify(var)
    b = certified_real_ball(f, x, domain)
    if b is not None:
        ok = {
            ">": b.lower.is_positive is True,
            "<": b.upper.is_negative is True,
            ">=": b.lower.is_nonnegative is True,
            "<=": b.upper.is_nonpositive is True,
        }.get(relation, False)
        if ok:
            return InequalityCertificate(
                f, x, domain, relation, TruthValue.TRUE, "outward-interval", b
            )
    # Taylor is used as a replayable tightening witness; unresolved contact with zero stays UNKNOWN.
    tm = taylor_model(f, x, domain, order=2)
    if max_depth > 0 and isinstance(domain, sp.Interval):
        mid = (domain.start + domain.end) / 2
        left = sp.Interval(domain.start, mid)
        right = sp.Interval(mid, domain.end)
        a = certify_inequality(
            lhs, rhs, x, left, relation=relation, max_depth=max_depth - 1
        )
        b2 = certify_inequality(
            lhs, rhs, x, right, relation=relation, max_depth=max_depth - 1
        )
        if a.classification is TruthValue.TRUE and b2.classification is TruthValue.TRUE:
            return InequalityCertificate(
                f,
                x,
                domain,
                relation,
                TruthValue.TRUE,
                "interval-subdivision",
                (a, b2, tm),
            )
    return InequalityCertificate(
        f, x, domain, relation, TruthValue.UNKNOWN, "inconclusive", tm
    )
