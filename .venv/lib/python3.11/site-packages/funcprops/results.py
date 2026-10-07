from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Protocol, runtime_checkable

from sympy import Basic, Piecewise, S, Symbol, sympify


class TruthValue(str, Enum):
    TRUE = "true"
    FALSE = "false"
    UNKNOWN = "unknown"


class Analyticity(str, Enum):
    ANALYTIC = "analytic"
    NOT_ANALYTIC = "not_analytic"
    UNKNOWN = "unknown"


class Meromorphicity(str, Enum):
    MEROMORPHIC = "meromorphic"
    NOT_MEROMORPHIC = "not_meromorphic"
    UNKNOWN = "unknown"


class Monotonicity(str, Enum):
    STRICTLY_INCREASING = "strictly_increasing"
    NONDECREASING = "nondecreasing"
    CONSTANT = "constant"
    NONINCREASING = "nonincreasing"
    STRICTLY_DECREASING = "strictly_decreasing"
    NOT_MONOTONE = "not_monotone"
    UNKNOWN = "unknown"


class Convexity(str, Enum):
    STRONGLY_CONVEX = "strongly_convex"
    STRICTLY_CONVEX = "strictly_convex"
    CONVEX = "convex"
    AFFINE = "affine"
    CONCAVE = "concave"
    STRICTLY_CONCAVE = "strictly_concave"
    STRONGLY_CONCAVE = "strongly_concave"
    LOG_CONVEX = "log_convex"
    LOG_CONCAVE = "log_concave"
    QUASICONVEX = "quasiconvex"
    QUASICONCAVE = "quasiconcave"
    PSEUDOCONVEX = "pseudoconvex"
    PSEUDOCONCAVE = "pseudoconcave"
    NEITHER = "neither"
    UNKNOWN = "unknown"


class Sign(str, Enum):
    POSITIVE = "positive"
    NONNEGATIVE = "nonnegative"
    ZERO = "zero"
    NONPOSITIVE = "nonpositive"
    NEGATIVE = "negative"
    SIGN_CHANGING = "sign_changing"
    UNKNOWN = "unknown"


class Boundedness(str, Enum):
    BOUNDED = "bounded"
    BOUNDED_ABOVE = "bounded_above"
    BOUNDED_BELOW = "bounded_below"
    UNBOUNDED = "unbounded"
    UNKNOWN = "unknown"


class Oscillation(str, Enum):
    NONOSCILLATORY = "nonoscillatory"
    FINITELY_OSCILLATORY = "finitely_oscillatory"
    OSCILLATORY = "oscillatory"
    UNKNOWN = "unknown"


class Periodicity(str, Enum):
    CONSTANT = "constant"
    PERIODIC = "periodic"
    APERIODIC = "aperiodic"
    UNKNOWN = "unknown"


class Symmetry(str, Enum):
    EVEN = "even"
    ODD = "odd"
    BOTH = "both"
    NEITHER = "neither"
    CENTRAL = "central"
    UNKNOWN = "unknown"


class Homogeneity(str, Enum):
    HOMOGENEOUS = "homogeneous"
    NOT_HOMOGENEOUS = "not_homogeneous"
    UNKNOWN = "unknown"


class Continuity(str, Enum):
    CONTINUOUS = "continuous"
    DISCONTINUOUS = "discontinuous"
    UNKNOWN = "unknown"


class Differentiability(str, Enum):
    SMOOTH = "smooth"
    CK = "Ck"
    DIFFERENTIABLE = "differentiable"
    NOT_DIFFERENTIABLE = "not_differentiable"
    UNKNOWN = "unknown"


class SingularityKind(str, Enum):
    POLE = "pole"
    ESSENTIAL = "essential"
    BRANCH_POINT = "branch_point"
    BRANCH_CUT = "branch_cut"
    DEFINITION_CUT = "definition_cut"
    REMOVABLE = "removable"
    DISCONTINUITY = "discontinuity"


@dataclass(frozen=True)
class ConditionalExpression:
    value: Any
    condition: Basic = S.true

    def __str__(self):
        return f"{self.condition} -> {self.value}"


@dataclass(frozen=True)
class Singularity:
    kind: SingularityKind
    condition: Basic
    guard: Basic = S.true
    multiplicity: int | Basic | None = None
    source: str | None = None

    @property
    def locus(self):
        return self.condition


@dataclass(frozen=True)
class SingularityResult:
    """Classified Boolean singularity loci with their supporting records.

    Every locus is a Boolean condition in the analysis variable(s).
    ``singularities`` is the classified union restricted to the effective function
    domain. ``unrestricted_singularities`` preserves the ambient singular locus,
    while ``interior_singularities``, ``boundary_singularities``, and
    ``exterior_singularities`` classify it relative to the restriction. ``records``
    preserves the lower-level :class:`Singularity` objects, including guards,
    multiplicities, and provenance.
    """

    singularities: Basic
    poles: Basic = S.false
    removable: Basic = S.false
    essential: Basic = S.false
    branch_points: Basic = S.false
    branch_cuts: Basic = S.false
    definition_cuts: Basic = S.false
    discontinuities: Basic = S.false
    records: tuple[Singularity, ...] = ()
    domain: Any = S.Complexes
    restriction: Basic = S.true
    unrestricted_singularities: Basic = S.false
    interior_singularities: Basic = S.false
    boundary_singularities: Basic = S.false
    exterior_singularities: Basic = S.false
    domain_interior: Basic = S.true
    domain_boundary: Basic = S.false
    domain_exterior: Basic = S.false

    @property
    def primary(self):
        return self.singularities

    @property
    def classification(self):
        return self.singularities

    @property
    def evidence(self):
        return self.records


@dataclass(frozen=True)
class PropertyResult:
    primary: Any
    cases: tuple[ConditionalExpression, ...] = ()
    all_properties: tuple[Any, ...] = ()
    source: str | None = None
    exceptions: tuple[Any, ...] = ()
    evidence: tuple[Any, ...] = ()

    @property
    def resolved(self):
        """Return the unconditional value when the cases imply one."""
        if self.primary is not None:
            return self.primary
        values = {case.value for case in self.cases if case.condition is S.true}
        return next(iter(values)) if len(values) == 1 else None

    def __str__(self):
        if self.primary is not None:
            return str(self.primary)
        if self.cases:
            return "; ".join(str(case) for case in self.cases)
        return "unknown"

    def _repr_html_(self):
        if self.primary is not None:
            return f"<b>Result:</b> {self.primary}"
        if self.cases:
            rows = "".join(
                f"<tr><td>{case.condition}</td><td>{case.value}</td></tr>"
                for case in self.cases
            )
            return (
                "<table><thead><tr><th>Condition</th><th>Result</th></tr></thead>"
                f"<tbody>{rows}</tbody></table>"
            )
        return "<b>Result:</b> unknown"

    def as_detailed_result(self):
        """Project this result onto the package-wide detailed-result protocol."""
        return as_detailed_result(self)

    def as_piecewise(self):
        """Return the compact symbolic value represented by this result.

        Conditional classifications become a SymPy ``Piecewise``.  An explicit
        ``unknown`` fallback is appended so an uncovered parameter region is
        never silently interpreted as false (or left to Piecewise's implicit
        ``nan`` fallback).  Exhaustive complementary cases are simplified by
        SymPy in the usual way.
        """
        if not self.cases:
            return self.primary

        def symbolic_value(value):
            if value is TruthValue.TRUE:
                return S.true
            if value is TruthValue.FALSE:
                return S.false
            if isinstance(value, Basic):
                return value
            if isinstance(value, Enum):
                return Symbol(str(value.value))
            return sympify(value, evaluate=False)

        pairs = [(symbolic_value(c.value), c.condition) for c in self.cases]
        pairs.append((Symbol("unknown"), S.true))
        return Piecewise(*pairs)


@dataclass(frozen=True)
class DetailedPropertyResult:
    """Uniform non-breaking view of detailed property-analysis output.

    Public property functions may keep returning their compact enums, sets, or
    specialized result dataclasses.  ``as_detailed_result`` projects any such
    value onto this common protocol for downstream orchestration.
    """

    primary: Any
    cases: tuple[ConditionalExpression, ...] = ()
    domain: Any = None
    complete: bool | None = None
    source: str | None = None
    evidence: tuple[Any, ...] = ()
    certificate: Any = None
    exceptions: tuple[Any, ...] = ()
    all_properties: tuple[Any, ...] = ()

    @property
    def resolved(self):
        """Return the compact value when this detailed result is unconditional."""
        if self.primary is not None:
            return self.primary
        values = {case.value for case in self.cases if case.condition is S.true}
        return next(iter(values)) if len(values) == 1 else None

    def __str__(self):
        return explain_result(self)

    def _repr_html_(self):
        fields = [
            ("Result", self.primary if self.primary is not None else "conditional")
        ]
        if self.domain is not None:
            fields.append(("Domain", self.domain))
        if self.cases:
            fields.append(("Cases", len(self.cases)))
        if self.source:
            fields.append(("Proof route", self.source.replace("-", " ")))
        fields.append(
            ("Certified", "yes" if self.certificate is not None else "not recorded")
        )
        rows = "".join(f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in fields)
        return f"<table><tbody>{rows}</tbody></table>"


@runtime_checkable
class SupportsDetailedPropertyResult(Protocol):
    """Structural protocol implemented by objects exposing detailed results."""

    def as_detailed_result(self) -> DetailedPropertyResult: ...


def as_detailed_result(value, *, domain=None, source=None, complete=None):
    """Return a uniform detailed view without changing compact public APIs."""
    if isinstance(value, DetailedPropertyResult):
        return value
    if isinstance(value, PropertyResult):
        certificate = next(
            (
                item
                for item in value.evidence
                if hasattr(item, "verify") or hasattr(item, "proof")
            ),
            None,
        )
        return DetailedPropertyResult(
            primary=value.primary,
            cases=value.cases,
            domain=domain,
            complete=complete,
            source=source or value.source,
            evidence=value.evidence,
            certificate=certificate,
            exceptions=value.exceptions,
            all_properties=value.all_properties,
        )
    primary = getattr(value, "classification", value)
    evidence = tuple(getattr(value, "evidence", ()))
    exceptions = tuple(getattr(value, "exceptional_loci", ()))
    result_domain = domain if domain is not None else getattr(value, "domain", None)
    return DetailedPropertyResult(
        primary=primary,
        cases=tuple(getattr(value, "cases", ())),
        domain=result_domain,
        complete=complete,
        source=source or getattr(value, "source", None),
        evidence=evidence,
        certificate=(
            getattr(value, "certificate", None)
            or (value if callable(getattr(value, "verify", None)) else None)
        ),
        exceptions=exceptions,
        all_properties=tuple(getattr(value, "all_properties", ())),
    )


def explain_result(value):
    """Return a concise explanation of a property result or certificate."""
    detail = as_detailed_result(value)
    parts = []
    if detail.primary is not None:
        parts.append(f"Result: {detail.primary}.")
    elif detail.cases:
        rendered = "; ".join(
            f"{case.condition} -> {case.value}" for case in detail.cases
        )
        parts.append(f"Cases: {rendered}.")
    else:
        parts.append("Result is unresolved.")
    if detail.source:
        parts.append(f"Proof route: {detail.source.replace('-', ' ')}.")
    if detail.domain is not None:
        parts.append(f"Domain: {detail.domain}.")
    if detail.certificate is not None:
        parts.append(f"Certificate: {type(detail.certificate).__name__}.")
    return " ".join(parts)
