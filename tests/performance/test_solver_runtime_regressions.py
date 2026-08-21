"""Opt-in runtime regression checks for representative CAD workloads.

Run with ``pytest -m performance``.  The budgets are deliberately generous so
these tests catch order-of-magnitude regressions rather than normal machine
variance.  Set ``SEMIALG_PERF_BUDGET_SCALE`` to scale all limits on slower or
faster CI workers.
"""

from __future__ import annotations

import os
import time

import pytest
import sympy as sp

from semialg import extract_cylindrical_solution, extract_explicit_cylindrical_solution

pytestmark = [pytest.mark.slow, pytest.mark.performance]


def _budget(seconds: float) -> float:
    scale = float(os.environ.get("SEMIALG_PERF_BUDGET_SCALE", "1"))
    if scale <= 0:
        raise ValueError("SEMIALG_PERF_BUDGET_SCALE must be positive")
    return seconds * scale


def test_explicit_three_dimensional_linear_cad_runtime():
    x, y, z = sp.symbols("x y z", real=True)
    formula = sp.And(x >= 0, x <= 1, y >= 0, y <= x, z >= 0, z <= x + y)

    start = time.perf_counter()
    solution = extract_explicit_cylindrical_solution(formula, (x, y, z))
    elapsed = time.perf_counter() - start

    assert solution is not None and solution.cells
    assert elapsed < _budget(8.0), f"explicit 3D CAD took {elapsed:.3f}s"


def test_nonlinear_two_dimensional_cad_runtime():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(x >= 0, x <= 1, y**2 <= x)

    start = time.perf_counter()
    solution = extract_cylindrical_solution(formula, (x, y))
    elapsed = time.perf_counter() - start

    assert solution.full_dimensional_cells
    assert elapsed < _budget(15.0), f"nonlinear 2D CAD took {elapsed:.3f}s"


def test_modular_boolean_oracle_runtime():
    from semialg.solve.integer.congruence import solve_quant_free_mod_sys

    x, y = sp.symbols("x y", integer=True)
    formula = sp.Or(sp.Ne(x + y, 0), sp.Eq(2 * x - y, 3))

    start = time.perf_counter()
    result = solve_quant_free_mod_sys(formula, (x, y), 30, max_points=2000)
    elapsed = time.perf_counter() - start

    assert result.complete and result.points
    assert elapsed < _budget(5.0), f"modular Boolean solve took {elapsed:.3f}s"


def test_periodic_transcendental_inequality_runtime():
    from semialg.solve.transcendental import build_trans_state, reduce_trans_problem

    x = sp.Symbol("x", real=True)
    state = build_trans_state(sp.sin(x) > 0, (x,))

    start = time.perf_counter()
    result = reduce_trans_problem(state)
    elapsed = time.perf_counter() - start

    assert result.method == "periodic_interval_reconstruction"
    assert result.result_semantics == "periodic_window_approximation"
    assert elapsed < _budget(5.0), f"periodic transcendental reduction took {elapsed:.3f}s"
