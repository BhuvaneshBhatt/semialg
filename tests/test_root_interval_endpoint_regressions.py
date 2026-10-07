"""Exact independent oracles for Sturm fallback endpoint conventions."""

from dataclasses import replace

import pytest
import sympy as sp

from semialg.algebraic import certify_polynomial_root_interval, verify_root_interval_certificate
from semialg.algebraic.roots import _count_roots_in_interval


@pytest.mark.parametrize("roots", [(-2, -1, 0), (-2, -1, 0, 1)])
@pytest.mark.parametrize("multiplicity", [1, 2])
@pytest.mark.parametrize(
    "include_left,include_right", [(False, False), (False, True), (True, False), (True, True)]
)
def test_sturm_certificate_counts_endpoints_once(roots, multiplicity, include_left, include_right):
    x = sp.Symbol("x")
    polynomial = sp.prod((x - r) ** multiplicity for r in roots)
    cert = certify_polynomial_root_interval(
        polynomial, -2, 1, var=x, include_left=include_left, include_right=include_right
    )
    expected = sum(
        -2 < r < 1 or (include_left and r == -2) or (include_right and r == 1) for r in roots
    )
    assert cert.method == "descartes+sturm"
    assert cert.root_count == expected
    assert verify_root_interval_certificate(cert)
    poly = sp.Poly(polynomial, x, domain=sp.QQ)
    assert (
        _count_roots_in_interval(
            poly,
            sp.Rational(-2),
            sp.Rational(1),
            exclude_left=not include_left,
            exclude_right=not include_right,
        )
        == expected
    )


@pytest.mark.parametrize("point", [-1, 2])
@pytest.mark.parametrize(
    "include_left,include_right", [(False, False), (False, True), (True, False), (True, True)]
)
def test_zero_width_sturm_count_never_double_subtracts(point, include_left, include_right):
    x = sp.Symbol("x")
    poly = sp.Poly((x + 1) ** 3, x, domain=sp.QQ)
    expected = int(point == -1 and include_left and include_right)
    assert (
        _count_roots_in_interval(
            poly,
            sp.Rational(point),
            sp.Rational(point),
            exclude_left=not include_left,
            exclude_right=not include_right,
        )
        == expected
    )
    cert = certify_polynomial_root_interval(
        poly, point, point, include_left=include_left, include_right=include_right
    )
    assert cert.root_count == expected
    assert verify_root_interval_certificate(cert)


def test_sturm_certificate_replay_rejects_incorrect_count():
    x = sp.Symbol("x")
    cert = certify_polynomial_root_interval(
        x * (x + 1) * (x + 2), -2, 1, var=x, include_left=False, include_right=False
    )
    assert cert.root_count == 2
    assert not verify_root_interval_certificate(replace(cert, root_count=1))
