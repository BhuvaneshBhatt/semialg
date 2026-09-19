"""Metamorphic properties for sparse SOS and CAD semantics."""

from __future__ import annotations

import pytest
import sympy as sp

from semialg.cad_algorithms.decomposition import decomp_collins_complete
from semialg.sos_certificates import sparse_sos_monomial_basis

pytestmark = pytest.mark.slow

hypothesis = pytest.importorskip("hypothesis")
given = hypothesis.given
settings = hypothesis.settings
st = pytest.importorskip("hypothesis.strategies")


@settings(max_examples=20, deadline=None)
@given(a=st.integers(1, 4), b=st.integers(1, 4), scale=st.integers(1, 4))
def test_sparse_sos_basis_is_invariant_under_positive_coefficient_scaling(a, b, scale):
    x, y = sp.symbols("x y")
    polynomial = a * x**8 + b * y**8
    assert sparse_sos_monomial_basis(polynomial, (x, y)) == sparse_sos_monomial_basis(
        scale * polynomial, (x, y)
    )


@settings(max_examples=16, deadline=None)
@given(a=st.integers(-3, 3), b=st.integers(-3, 3))
def test_cad_truth_is_invariant_under_factorized_or_expanded_input(a, b):
    x, y = sp.symbols("x y")
    source = (x - a) * (x - b) * (y - 1)
    factored = decomp_collins_complete((source,), (x, y))
    expanded = decomp_collins_complete((sp.expand(source),), (x, y))
    assert factored.complete and expanded.complete
    assert factored.cell_count_by_level() == expanded.cell_count_by_level()
