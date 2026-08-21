from __future__ import annotations


class SemialgError(Exception):
    """Base class for semialg-specific failures."""


class SemialgStrategyFailure(SemialgError):
    """A speculative exact/symbolic strategy could not handle the input."""


class ReconstructionFailure(SemialgStrategyFailure):
    """A symbolic reconstruction strategy could not certify a representation."""


class DimensionMismatchError(SemialgError, ValueError):
    """Parallel coordinates, variables, or other dimensions do not agree."""


class CertificationFailure(SemialgError):
    """A requested mathematical certificate could not be established."""


class ExactEvaluationFailure(SemialgStrategyFailure):
    """An exact backend could not evaluate an otherwise valid object."""


__all__ = [
    "SemialgError",
    "SemialgStrategyFailure",
    "ReconstructionFailure",
    "DimensionMismatchError",
    "CertificationFailure",
    "ExactEvaluationFailure",
]
