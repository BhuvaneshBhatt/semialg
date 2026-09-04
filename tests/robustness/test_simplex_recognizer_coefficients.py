from __future__ import annotations

import sympy as sp

from semialg.region_integrate import _reduce_standard_simplex_2d


def test_simplex_recognizer_uses_canonical_linear_coefficients() -> None:
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x >= 0, y >= 0, x + y <= 1)

    pieces = _reduce_standard_simplex_2d(x * y, condition, x, y)

    assert pieces is not None
    assert len(pieces) == 1
    assert pieces[0].method == "unit_simplex_iterated_integral"


def test_simplex_recognizer_accepts_reversed_equivalent_relations() -> None:
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(-x <= 0, -y <= 0, 1 - x - y >= 0)

    pieces = _reduce_standard_simplex_2d(1, condition, x, y)

    assert pieces is not None
    assert pieces[0].method == "unit_simplex_iterated_integral"
