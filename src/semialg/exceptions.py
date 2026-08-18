"""Shared exception hierarchy for semialgebraic algorithms."""

from __future__ import annotations


class SemialgError(Exception):
    """Base class for package-specific failures."""


class UnsupportedFragmentError(SemialgError, ValueError):
    """Raised when an input is outside a supported symbolic fragment."""


class BackendFailure(SemialgError):
    """Raised when an optional backend fails after accepting an input."""


class FormulaNormalizationError(SemialgError, ValueError):
    """Raised when a formula cannot be normalized for a solver."""


class AlgebraicSolvingError(SemialgError, ValueError):
    """Raised by exact algebraic solving backends."""


class QuantifierEliminationError(SemialgError, ValueError):
    """Raised by quantifier-elimination backends."""
