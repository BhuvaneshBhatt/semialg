"""Diagnostic helpers for semialgebraic solving results."""

from __future__ import annotations

import sympy as sp


def solution_capability_diagnostics(input_formula: sp.Expr, **extra: object) -> dict[str, object]:
    """Return stable capability-oriented diagnostics for solve results."""

    diagnostics: dict[str, object] = {
        "input_formula": sp.sstr(input_formula),
        "solver_stage": "semialgebraic_solution",
        "capabilities": {
            "component_decomposition": "one_dimensional",
            "parameter_conditions": True,
            "vertical_cells": "two_dimensional",
            "structural_sampling": True,
            "set_operations": True,
            "plotting_data": True,
        },
        "full_cylindrical_solution": True,
    }
    diagnostics.update(extra)
    return diagnostics
