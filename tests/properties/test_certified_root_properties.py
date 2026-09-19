"""Property-based contracts for exact certified polynomial-root utilities."""

from __future__ import annotations

import pytest
import sympy as sp
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg.algebraic import (
    certify_polynomial_root_interval,
    isolate_real_roots,
    rational_between_algebraic_reals,
    verify_polynomial_root_interval_certificate,
)

pytestmark = pytest.mark.slow


@settings(max_examples=40, deadline=None)
@given(
    roots=st.lists(st.integers(-5, 5), min_size=1, max_size=5),
    left=st.integers(-6, 4),
    width=st.integers(0, 6),
    include_left=st.booleans(),
    include_right=st.booleans(),
)
def test_interval_certificate_counts_distinct_generated_rational_roots(
    roots, left, width, include_left, include_right
):
    x = sp.Symbol("x")
    right = left + width
    polynomial = sp.prod(x - root for root in roots)
    certificate = certify_polynomial_root_interval(
        polynomial,
        left,
        right,
        var=x,
        include_left=include_left,
        include_right=include_right,
    )

    expected = {
        root
        for root in roots
        if (left < root < right)
        or (include_left and root == left)
        or (include_right and root == right)
    }
    if left == right and not (include_left and include_right):
        expected = set()

    assert certificate.root_count == len(expected)
    assert verify_polynomial_root_interval_certificate(certificate)


@settings(max_examples=20, deadline=None)
@given(n=st.sampled_from((2, 3, 5, 6, 7, 10, 11, 13)))
def test_separator_is_exactly_between_quadratic_conjugate_roots(n):
    x = sp.Symbol("x")
    left, right = isolate_real_roots(sp.Poly(x**2 - n, x, domain=sp.QQ))
    separator = rational_between_algebraic_reals(left, right)

    assert separator.is_Rational
    assert sp.simplify(left.as_expr() < separator) is sp.true
    assert sp.simplify(separator < right.as_expr()) is sp.true


@settings(max_examples=25, deadline=None)
@given(root=st.integers(-4, 4), multiplicity=st.integers(1, 4))
def test_root_interval_count_is_invariant_under_repeated_factors(root, multiplicity):
    x = sp.Symbol("x")
    polynomial = (x - root) ** multiplicity
    certificate = certify_polynomial_root_interval(polynomial, root - 1, root + 1, var=x)

    assert certificate.root_count == 1
    assert certificate.unique
