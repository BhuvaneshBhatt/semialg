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


def test_simplex_recognizer_is_invariant_under_positive_scaling() -> None:
    x, y = sp.symbols("x y", real=True)
    canonical = sp.And(x >= 0, y >= 0, x + y <= 1)
    scaled = sp.And(2 * x >= 0, 3 * y >= 0, 5 * x + 5 * y <= 5)

    canonical_pieces = _reduce_standard_simplex_2d(x * y, canonical, x, y)
    scaled_pieces = _reduce_standard_simplex_2d(x * y, scaled, x, y)

    assert canonical_pieces is not None and scaled_pieces is not None
    assert canonical_pieces[0].limits == scaled_pieces[0].limits
    assert canonical_pieces[0].method == scaled_pieces[0].method
