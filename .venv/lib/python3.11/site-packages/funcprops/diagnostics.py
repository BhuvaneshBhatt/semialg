"""Uniform unresolved-proof and completeness vocabulary."""

from enum import Enum


class UnknownReason(str, Enum):
    UNRESOLVED_REALNESS = "unresolved-realness"
    UNRESOLVED_SIGN = "unresolved-sign"
    INCOMPLETE_ZERO_GEOMETRY = "incomplete-zero-geometry"
    UNSUPPORTED_TRANSCENDENTAL = "unsupported-transcendental"
    PARAMETER_STRATIFICATION_INCOMPLETE = "parameter-stratification-incomplete"
    BACKEND_INCONCLUSIVE = "backend-inconclusive"
    PROOF_BUDGET_EXHAUSTED = "proof-budget-exhausted"
    DOMAIN_GEOMETRY_UNSUPPORTED = "domain-geometry-unsupported"


class Completeness(str, Enum):
    COMPLETE = "complete"
    COMPLETE_ON_DOMAIN = "complete-on-domain"
    LOCAL = "local"
    PARTIAL = "partial"
    LOWER_BOUND = "lower-bound"
    UNKNOWN = "unknown"
