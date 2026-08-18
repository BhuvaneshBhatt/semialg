"""Shared status values for solver result objects."""

from __future__ import annotations

from enum import Enum


class SolverStatus(str, Enum):
    """Common satisfiability/status vocabulary."""

    SAT = "sat"
    UNSAT = "unsat"
    UNKNOWN = "unknown"
    PARTIAL = "partial"
    ERROR = "error"


class CoverageStatus(str, Enum):
    """Whether a specialized backend covered the requested problem."""

    COMPLETE = "complete"
    PARTIAL = "partial"
    NONE = "none"
