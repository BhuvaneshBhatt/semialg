"""Structural function geometry, boundary behavior, and attainment certificates.

``FunctionGeometry`` combines regularity, convexity, singularity, boundary escape,
level-boundedness, and structural mode-existence information. These routines are
proof-conservative: unsupported geometry is reported as ``TruthValue.UNKNOWN``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp
from sympy import S

from .exceptions import SYMBOLIC_QUERY_ERRORS
from .real_analysis.common import normalize_domain_set, normalize_variables
from .real_analysis.convexity import convexity as _convexity
from .real_analysis.smoothness import smoothness as _smoothness
from .results import Analyticity, Continuity, Convexity, SingularityResult, TruthValue
from .singularities import (
    analytic as _analytic,
)
from .singularities import (
    continuous,
)
from .singularities import (
    singularities as _singularities,
)


def _real_set(domain, variable):
    return normalize_domain_set(domain, variable)


def _is_bounded_real_set(domain):
    if isinstance(domain, sp.Interval):
        return domain.start is not S.NegativeInfinity and domain.end is not S.Infinity
    if isinstance(domain, sp.FiniteSet):
        return True
    if isinstance(domain, sp.Union):
        return all(_is_bounded_real_set(piece) for piece in domain.args)
    bounded = getattr(domain, "is_bounded", None)
    return bounded is True


def _set_closed(domain):
    if domain is S.Reals:
        return TruthValue.TRUE
    if isinstance(domain, sp.Set):
        closed = getattr(domain, "is_closed", None)
        if closed is True:
            return TruthValue.TRUE
        if closed is False:
            return TruthValue.FALSE
        if isinstance(domain, sp.Interval):
            return (
                TruthValue.TRUE
                if not domain.left_open and not domain.right_open
                else TruthValue.FALSE
            )
    return TruthValue.UNKNOWN


@dataclass(frozen=True)
class BoundaryLimit:
    endpoint: Any
    value: Any
    direction: str
    certified: bool
    source: str = "symbolic-limit"


@dataclass(frozen=True)
class BoundaryBehavior:
    domain: Any
    limits: tuple[BoundaryLimit, ...]
    complete: bool
    source: str

    @property
    def positive_infinity_on_escape(self):
        if not self.complete:
            return TruthValue.UNKNOWN
        escape = [
            item
            for item in self.limits
            if item.endpoint in (S.NegativeInfinity, S.Infinity)
        ]
        if not escape:
            return TruthValue.TRUE
        if all(item.certified and item.value is S.Infinity for item in escape):
            return TruthValue.TRUE
        if any(item.certified and item.value is not S.Infinity for item in escape):
            return TruthValue.FALSE
        return TruthValue.UNKNOWN


@dataclass(frozen=True)
class LevelBoundednessCertificate:
    classification: TruthValue
    domain: Any
    boundary_behavior: BoundaryBehavior | None
    source: str
    evidence: tuple[Any, ...] = ()

    @property
    def primary(self):
        return self.classification

    def verify(self):
        if self.classification is TruthValue.TRUE:
            if self.source == "bounded-domain":
                return (
                    isinstance(self.domain, sp.Set)
                    and self.domain is not S.Reals
                    and _is_bounded_real_set(self.domain)
                )
            if self.source == "boundary-plus-infinity":
                return (
                    self.boundary_behavior is not None
                    and self.boundary_behavior.positive_infinity_on_escape
                    is TruthValue.TRUE
                )
            if self.source == "strong-convexity":
                return any(item is Convexity.STRONGLY_CONVEX for item in self.evidence)
            if self.source == "separable-coercivity":
                return bool(self.evidence) and all(
                    isinstance(item, LevelBoundednessCertificate)
                    and item.classification is TruthValue.TRUE
                    and item.verify()
                    for item in self.evidence
                )
            return False
        if self.classification is TruthValue.FALSE:
            if (
                self.source != "escaping-nonpositive-infinity"
                or self.boundary_behavior is None
            ):
                return False
            return any(
                item.certified
                and item.endpoint in (S.NegativeInfinity, S.Infinity)
                and item.value is not S.Infinity
                for item in self.boundary_behavior.limits
            )
        return self.source in {"insufficient-structure", "multivariate-unresolved"}


@dataclass(frozen=True)
class ModeExistenceCertificate:
    sense: str
    exists: bool | None
    domain: Any
    continuity: Continuity
    level_boundedness: LevelBoundednessCertificate
    domain_closed: TruthValue
    source: str
    evidence: tuple[Any, ...] = ()

    @property
    def classification(self):
        if self.exists is True:
            return TruthValue.TRUE
        if self.exists is False:
            return TruthValue.FALSE
        return TruthValue.UNKNOWN

    @property
    def primary(self):
        return self.classification

    def verify(self):
        if self.exists is True:
            return (
                self.source == "continuous-level-bounded-weierstrass"
                and self.continuity is Continuity.CONTINUOUS
                and self.domain_closed is TruthValue.TRUE
                and self.level_boundedness.classification is TruthValue.TRUE
                and self.level_boundedness.verify()
            )
        if self.exists is None:
            return self.source == "insufficient-structure"
        return False


@dataclass(frozen=True)
class FunctionGeometry:
    expression: Any
    variables: tuple[Any, ...]
    domain: Any
    continuity: Continuity
    ck_order: Any
    c2: TruthValue
    analytic: Analyticity
    convexity: Convexity
    singularities: SingularityResult
    boundary_behavior: BoundaryBehavior
    level_boundedness: LevelBoundednessCertificate
    minimum_exists: ModeExistenceCertificate
    maximum_exists: ModeExistenceCertificate
    evidence: tuple[Any, ...] = ()


def _unbounded_endpoints(domain):
    if domain is S.Reals:
        return ((S.NegativeInfinity, "+"), (S.Infinity, "-"))
    if isinstance(domain, sp.Interval):
        out = []
        if domain.start is S.NegativeInfinity:
            out.append((S.NegativeInfinity, "+"))
        if domain.end is S.Infinity:
            out.append((S.Infinity, "-"))
        return tuple(out)
    if isinstance(domain, sp.Union):
        out = []
        for piece in domain.args:
            out.extend(_unbounded_endpoints(piece))
        # deterministic de-duplication
        return tuple(dict.fromkeys(out))
    return ()


def boundary_behavior(expr, var, *, domain=S.Reals):
    """Certify explicit univariate limits at unbounded real-domain endpoints."""
    expr = sp.sympify(expr)
    vv = normalize_variables(var)
    if len(vv) != 1:
        return BoundaryBehavior(domain, (), False, "multivariate-radial-unresolved")
    variable = vv[0]
    real_domain = _real_set(domain, variable)
    if real_domain is None:
        return BoundaryBehavior(domain, (), False, "domain-unresolved")
    endpoints = _unbounded_endpoints(real_domain)
    if not endpoints:
        complete = _is_bounded_real_set(real_domain)
        return BoundaryBehavior(
            real_domain,
            (),
            complete,
            "bounded-domain" if complete else "no-unbounded-endpoints",
        )
    limits = []
    complete = True
    for endpoint, direction in endpoints:
        try:
            value = sp.limit(expr, variable, endpoint, dir=direction)
            certified = not isinstance(value, sp.Limit)
        except SYMBOLIC_QUERY_ERRORS:
            value = sp.Limit(expr, variable, endpoint, dir=direction)
            certified = False
        limits.append(BoundaryLimit(endpoint, value, direction, certified))
        complete = complete and certified
    return BoundaryBehavior(real_domain, tuple(limits), complete, "symbolic-limits")


def _separable_terms(expr, variables):
    terms = sp.Add.make_args(sp.expand(expr))
    buckets = {v: [] for v in variables}
    constant = []
    for term in terms:
        deps = [v for v in variables if term.has(v)]
        if not deps:
            constant.append(term)
        elif len(deps) == 1:
            buckets[deps[0]].append(term)
        else:
            return None
    return tuple(sp.Add(*buckets[v]) for v in variables)


def _univariate_level_boundedness(expr, variable, domain, behavior=None):
    real_domain = _real_set(domain, variable)
    if behavior is None:
        behavior = boundary_behavior(expr, variable, domain=domain)
    if (
        isinstance(real_domain, sp.Set)
        and real_domain is not S.Reals
        and _is_bounded_real_set(real_domain)
    ):
        return LevelBoundednessCertificate(
            TruthValue.TRUE, real_domain, behavior, "bounded-domain", (real_domain,)
        )
    verdict = behavior.positive_infinity_on_escape
    if verdict is TruthValue.TRUE:
        return LevelBoundednessCertificate(
            TruthValue.TRUE,
            real_domain,
            behavior,
            "boundary-plus-infinity",
            behavior.limits,
        )
    if verdict is TruthValue.FALSE and behavior.complete:
        return LevelBoundednessCertificate(
            TruthValue.FALSE,
            real_domain,
            behavior,
            "escaping-nonpositive-infinity",
            behavior.limits,
        )
    return LevelBoundednessCertificate(
        TruthValue.UNKNOWN, real_domain, behavior, "insufficient-structure"
    )


def _resolved_domain_and_closed(domain, variables):
    if len(variables) == 1:
        resolved = _real_set(domain, variables[0])
        return resolved, _set_closed(resolved)
    if domain in (S.Reals, True, S.true, None):
        return S.Reals, TruthValue.TRUE
    return domain, TruthValue.UNKNOWN


def _mode_existence_from_facts(
    sense, continuity, level_boundedness, resolved_domain, domain_closed
):
    if (
        continuity is Continuity.CONTINUOUS
        and level_boundedness.classification is TruthValue.TRUE
        and domain_closed is TruthValue.TRUE
    ):
        return ModeExistenceCertificate(
            sense,
            True,
            resolved_domain,
            continuity,
            level_boundedness,
            domain_closed,
            "continuous-level-bounded-weierstrass",
            (level_boundedness,),
        )
    return ModeExistenceCertificate(
        sense,
        None,
        resolved_domain,
        continuity,
        level_boundedness,
        domain_closed,
        "insufficient-structure",
    )


def _level_boundedness_certificate(expr, variables, *, domain=S.Reals):
    """Return whether all sublevel sets are bounded."""
    expr = sp.sympify(expr)
    vv = normalize_variables(variables)
    if len(vv) == 1:
        return _univariate_level_boundedness(expr, vv[0], domain)

    # Current multivariate certificates target the whole Euclidean domain.
    whole_space = domain in (S.Reals, True, S.true, None)
    if whole_space:
        convexity = _convexity(expr, vv)
        if convexity is Convexity.STRONGLY_CONVEX:
            return LevelBoundednessCertificate(
                TruthValue.TRUE, S.Reals, None, "strong-convexity", (convexity,)
            )
        separable = _separable_terms(expr, vv)
        if separable is not None:
            component_certs = tuple(
                _level_boundedness_certificate(term, var)
                for term, var in zip(separable, vv, strict=True)
            )
            if all(cert.classification is TruthValue.TRUE for cert in component_certs):
                return LevelBoundednessCertificate(
                    TruthValue.TRUE,
                    S.Reals,
                    None,
                    "separable-coercivity",
                    component_certs,
                )
    return LevelBoundednessCertificate(
        TruthValue.UNKNOWN, domain, None, "multivariate-unresolved"
    )


def level_boundedness(expr, variables, *, domain=S.Reals, return_result=False):
    """Return whether all sublevel sets are bounded.

    Use ``return_result=True`` for the structural certificate and proof source.
    """
    certificate = _level_boundedness_certificate(expr, variables, domain=domain)
    return certificate if return_result else certificate.classification


def _mode_existence_certificate(expr, variables, *, domain=S.Reals, sense="max"):
    """Certify existence of a global minimum or maximum from structural facts."""
    if sense not in {"min", "max"}:
        raise ValueError("sense must be 'min' or 'max'")
    expr = sp.sympify(expr)
    vv = normalize_variables(variables)
    objective = expr if sense == "min" else -expr
    level_boundedness = _level_boundedness_certificate(objective, vv, domain=domain)
    continuity = continuous(expr, vv if len(vv) > 1 else vv[0], domain=domain)
    resolved_domain, domain_closed = _resolved_domain_and_closed(domain, vv)
    return _mode_existence_from_facts(
        sense, continuity, level_boundedness, resolved_domain, domain_closed
    )


def mode_existence(
    expr, variables, *, domain=S.Reals, sense="max", return_result=False
):
    """Return whether a global minimum or maximum is structurally guaranteed.

    Use ``return_result=True`` for the attainment certificate and supporting facts.
    """
    certificate = _mode_existence_certificate(
        expr, variables, domain=domain, sense=sense
    )
    return certificate if return_result else certificate.classification


def geometry(expr, variables, *, domain=S.Reals, assumptions=S.true):
    """Summarize regularity, convexity, singularity, and escape geometry."""
    expr = sp.sympify(expr)
    vv = normalize_variables(variables)
    var_arg = vv if len(vv) > 1 else vv[0]
    continuity = continuous(expr, var_arg, domain=domain, assumptions=assumptions)
    analytic = _analytic(expr, var_arg, domain=domain, assumptions=assumptions)
    if len(vv) == 1:
        smooth = _smoothness(
            expr, vv[0], domain=domain, max_order=2, return_result=True
        )
        ck_order = getattr(smooth, "highest_continuous_derivative", None)
        if (
            analytic is Analyticity.ANALYTIC
            or ck_order is S.Infinity
            or ck_order == sp.oo
        ):
            c2 = TruthValue.TRUE
        elif isinstance(ck_order, int) and ck_order >= 2:
            c2 = TruthValue.TRUE
        else:
            c2 = TruthValue.UNKNOWN
    else:
        ck_order = sp.oo if analytic is Analyticity.ANALYTIC else None
        c2 = TruthValue.TRUE if analytic is Analyticity.ANALYTIC else TruthValue.UNKNOWN
    convexity = _convexity(expr, var_arg, domain=domain)
    singularities = _singularities(
        expr, var_arg, domain=S.Complexes, assumptions=assumptions, return_result=True
    )
    boundary = boundary_behavior(expr, var_arg, domain=domain)
    if len(vv) == 1:
        level = _univariate_level_boundedness(expr, vv[0], domain, boundary)
    else:
        level = _level_boundedness_certificate(expr, vv, domain=domain)
    resolved_domain, domain_closed = _resolved_domain_and_closed(domain, vv)
    minimum = _mode_existence_from_facts(
        "min", continuity, level, resolved_domain, domain_closed
    )
    maximum_level = _level_boundedness_certificate(-expr, vv, domain=domain)
    maximum = _mode_existence_from_facts(
        "max", continuity, maximum_level, resolved_domain, domain_closed
    )
    return FunctionGeometry(
        expr,
        vv,
        domain,
        continuity,
        ck_order,
        c2,
        analytic,
        convexity,
        singularities,
        boundary,
        level,
        minimum,
        maximum,
        (singularities, level),
    )
