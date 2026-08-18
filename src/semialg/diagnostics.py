from __future__ import annotations

from .invariants import validate_cad_result
from .model import CADResult


def summarize_diagnostics(cad: CADResult) -> dict:
    invariant_issues = validate_cad_result(cad)
    return {
        "projection_mode": cad.projection_mode,
        "well_oriented": cad.well_oriented,
        "used_fallback_projection": cad.used_fallback_projection,
        "nullification_events": len(cad.nullification_events),
        "formula_pruned_cells": len(cad.formula_pruned_cells),
        "events": len(cad.diagnostics.events),
        "counters": dict(cad.diagnostics.counters),
        "cache_stats": {
            name: {"hits": stats.hits, "misses": stats.misses}
            for name, stats in cad.diagnostics.cache_stats.items()
        },
        "invariant_issues": invariant_issues,
    }
