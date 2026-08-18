"""Compatibility-free public equivalence helpers.

The validation package is the canonical home for equivalence and
symmetric-difference checks; this top-level module is a small public facade.
"""

from __future__ import annotations

from .validation.equivalence import EquivalenceReport, EquivCounterex, sym_diff_empty

__all__ = ["EquivCounterex", "EquivalenceReport", "sym_diff_empty"]
