"""Numerical optimizer hints become evidence only after exact Gram reconstruction."""

import pytest
import sympy as sp

from semialg.sdp_backends import external_sdp_search, reconstruct_exact_gram
from semialg.sos_certificates import SOSCertificate, verify_sos_certificate


def test_optimizer_hint():
    cvxpy = pytest.importorskip("cvxpy")
    if "CLARABEL" not in cvxpy.installed_solvers():
        pytest.skip("Clarabel is not installed")
    x, y = sp.symbols("x y")
    polynomial = x**4 + y**4
    candidate = external_sdp_search(polynomial, (x, y), backend="clarabel")
    assert candidate.success
    gram = reconstruct_exact_gram(
        polynomial, (x, y), candidate.monomial_basis, candidate.gram_matrix, max_denominator=1000
    )
    assert gram is not None
    certificate = SOSCertificate(polynomial, (x, y), candidate.monomial_basis, gram)
    assert verify_sos_certificate(polynomial, certificate, (x, y))
    vector = sp.Matrix(candidate.monomial_basis)
    assert sp.expand((vector.T * gram * vector)[0] - polynomial) == 0
