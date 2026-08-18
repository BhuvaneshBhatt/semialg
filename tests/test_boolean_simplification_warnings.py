from __future__ import annotations

import warnings

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

from semialg import solve_semialgebraic


def test_triangle_solve_emits_no_sympy_boolean_mul_deprecation_warning():
    x, y = sp.symbols("x y", real=True)
    triangle = sp.And(x >= 0, x <= 1, y >= x, y <= 1)

    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always", SymPyDeprecationWarning)
        solution = solve_semialgebraic(triangle, [x, y], count=0)

    assert solution.satisfiable is True
    assert not [item for item in recorded if issubclass(item.category, SymPyDeprecationWarning)]
