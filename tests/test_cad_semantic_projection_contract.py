"""Representation-independent contracts for factorized CAD projection."""

from __future__ import annotations

import sympy as sp

from semialg.cad_algorithms.decomposition import decomp_collins_complete
from tests._support.cad_assertions import sign_from_factor_signs


def test_source_polynomial_sign_is_recoverable_from_factorized_sign_table():
    x, y = sp.symbols("x y")
    source = (x - 1) * (x + 1) * (y - 1)
    cad = decomp_collins_complete((source,), (x, y))

    for cell in cad.cells:
        reconstructed = sign_from_factor_signs(source, (x, y), cell.signs)
        sample = dict(zip((x, y), cell.sample_exprs, strict=True))
        actual = int(sp.sign(sp.expand(source).subs(sample)))
        assert reconstructed == actual
