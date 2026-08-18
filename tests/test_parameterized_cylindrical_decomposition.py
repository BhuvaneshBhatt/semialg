from __future__ import annotations

import pytest
import sympy as sp

from semialg import parameterized_cylindrical_decomposition, solve_semialgebraic

pytestmark = pytest.mark.slow


def test_parameterized_cylindrical_decomposition_for_quadratic_feasibility():
    x, a, b = sp.symbols("x a b", real=True)
    decomp = parameterized_cylindrical_decomposition([sp.Eq(x**2 + a * x + b, 0)], [x], [a, b])
    assert decomp.parameter_condition != sp.false
    assert decomp.strata
    assert all(stratum.condition != sp.false for stratum in decomp.strata)


def test_solve_semialgebraic_exposes_parameter_decomposition_output():
    x, a = sp.symbols("x a", real=True)
    decomp = solve_semialgebraic(
        [x**2 + a < 0], [x], parameters=[a], output="parameter_decomposition"
    )
    assert decomp is not None
    assert sp.simplify_logic(decomp.parameter_condition) == (a < 0)
    assert decomp.stratum_count >= 1


def test_solve_semialgebraic_result_records_parameter_decomposition():
    x, a = sp.symbols("x a", real=True)
    result = solve_semialgebraic([x**2 + a < 0], [x], parameters=[a], count=0)
    assert result.parameter_conditions == (a < 0)
    assert result.parameter_decomposition is not None
