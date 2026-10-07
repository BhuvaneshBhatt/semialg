"""Certified complete monotonicity on real intervals."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..results import ConditionalExpression, PropertyResult, TruthValue
from ..statistical_special_functions import MittagLefflerE
from ..stratification import normalize_cases, semantic_condition
from .common import condition_to_set


@dataclass(frozen=True)
class CompleteMonotonicityCertificate:
    expression: sp.Basic
    variable: sp.Symbol
    domain: object
    classification: TruthValue
    assumptions: sp.Basic = S.true
    source: str = "insufficient-evidence"
    evidence: tuple = ()

    @property
    def primary(self):
        return self.classification

    def verify(self):
        value = complete_monotonicity(
            self.expression,
            self.variable,
            domain=self.domain,
            assumptions=self.assumptions,
            return_result=False,
        )
        if isinstance(value, PropertyResult):
            return any(
                c.value is self.classification
                and sp.simplify_logic(sp.Implies(self.assumptions, c.condition))
                is S.true
                for c in value.cases
            )
        return value is self.classification


def _nonnegative_axis(domain, var):
    ds = condition_to_set(domain, var)
    if isinstance(ds, sp.Interval) and ds.start == 0 and ds.end is S.Infinity:
        return True
    try:
        return (
            sp.simplify_logic(sp.Implies(semantic_condition(domain), var >= 0))
            is S.true
        )
    except (TypeError, ValueError, NotImplementedError):
        return False


def _positive_axis(domain, var):
    ds = condition_to_set(domain, var)
    if (
        isinstance(ds, sp.Interval)
        and ds.start == 0
        and ds.left_open
        and ds.end is S.Infinity
    ):
        return True
    try:
        return (
            sp.simplify_logic(sp.Implies(semantic_condition(domain), var > 0)) is S.true
        )
    except (TypeError, ValueError, NotImplementedError):
        return False


def complete_monotonicity(
    expr, var, *, domain=S.Reals, assumptions=S.true, return_result=False
):
    """Certify ``(-1)^n f^(n)(x) >= 0`` for every integer ``n>=0``.

    The implementation is theorem-driven: it returns UNKNOWN
    finite number of derivatives and extrapolating.
    """
    expr, var = sp.sympify(expr), sp.sympify(var)
    assumptions = sp.sympify(assumptions)
    classification = TruthValue.UNKNOWN
    source = "insufficient-evidence"
    cases = ()

    if _nonnegative_axis(domain, var):
        # Pollard's theorem: E_alpha(-x) is completely monotone for 0<alpha<=1.
        if expr.func is MittagLefflerE and len(expr.args) == 2 and expr.args[1] == -var:
            alpha = expr.args[0]
            guard = sp.And(alpha > 0, alpha <= 1)
            cases = (ConditionalExpression(TruthValue.TRUE, guard),)
            source = "pollard-mittag-leffler"
        elif expr.func is sp.exp and expr.args[0] == -var:
            classification = TruthValue.TRUE
            source = "elementary-exponential"
    if _positive_axis(domain, var):
        if expr.func is sp.besselk and len(expr.args) == 2 and expr.args[1] == var:
            nu = expr.args[0]
            cases = (ConditionalExpression(TruthValue.TRUE, sp.Contains(nu, S.Reals)),)
            source = "besselk-integral-representation"

    if cases:
        out = normalize_cases(cases, source=source)
        if out.primary is not None:
            classification = out.primary
        elif assumptions is not S.true:
            for c in out.cases:
                try:
                    if (
                        sp.simplify_logic(sp.Implies(assumptions, c.condition))
                        is S.true
                    ):
                        classification = c.value
                        break
                except (TypeError, ValueError, NotImplementedError):
                    pass
        if classification is TruthValue.UNKNOWN and not return_result:
            return out

    cert = CompleteMonotonicityCertificate(
        expr, var, domain, classification, assumptions, source, cases
    )
    return cert if return_result else classification
