"""Public exception exports for semialg.

The canonical exception classes are defined in :mod:`semialg.errors` and are
re-exported here as a compact import surface.
"""

from __future__ import annotations

from .errors import (
    AlgebraicSolvingError,
    BackendFailure,
    CertificationFailure,
    DimensionMismatchError,
    ExactEvaluationFailure,
    FormulaNormalizationError,
    QuantifierEliminationError,
    ReconstructionFailure,
    ResourceLimitError,
    SemialgError,
    SemialgStrategyFailure,
    UnsupportedFragmentError,
)

__all__ = [
    "SemialgError",
    "SemialgStrategyFailure",
    "UnsupportedFragmentError",
    "BackendFailure",
    "FormulaNormalizationError",
    "AlgebraicSolvingError",
    "QuantifierEliminationError",
    "ReconstructionFailure",
    "DimensionMismatchError",
    "CertificationFailure",
    "ExactEvaluationFailure",
    "ResourceLimitError",
]
