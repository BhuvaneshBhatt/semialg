import pytest
import sympy as sp

from semialg.constraint_certificates import (
    NonnegativeCombinationCertificate,
    verify_nonnegative_certificate,
)
from semialg.sos_certificates import SOSCertificate, verify_psd_exact, verify_sos_certificate

x = sp.Symbol("x", real=True)


@pytest.mark.parametrize(
    "polynomial,basis", [(-1, (sp.I,)), (-(x**2), (sp.I * x,)), (x, (sp.sqrt(x),))]
)
def test_invalid_sos_basis(polynomial, basis):
    certificate = SOSCertificate(sp.sympify(polynomial), (x,), basis, sp.ImmutableMatrix([[1]]))
    assert not verify_sos_certificate(polynomial, certificate)
    combination = NonnegativeCombinationCertificate(
        sp.sympify(polynomial), (), (), sp.sympify(polynomial), certificate
    )
    assert not verify_nonnegative_certificate(combination, (x,))


def test_complex_symmetric_matrix_is_not_psd():
    matrix = sp.ImmutableMatrix([[1, sp.I], [sp.I, -sp.Rational(1, 2)]])
    assert not verify_psd_exact(matrix).verified


def test_real_algebraic_polynomial_basis():
    polynomial = (sp.sqrt(2) * x + 1) ** 2
    certificate = SOSCertificate(polynomial, (x,), (sp.sqrt(2) * x + 1,), sp.ImmutableMatrix([[1]]))
    assert verify_sos_certificate(polynomial, certificate)


@pytest.mark.parametrize("basis,gram", [((sp.Float(1) * x,), [[1]]), ((x,), [[sp.Float(1)]])])
def test_inexact_proof_data(basis, gram):
    certificate = SOSCertificate(x**2, (x,), basis, sp.ImmutableMatrix(gram))
    assert not verify_sos_certificate(x**2, certificate)


def test_constant_sos():
    certificate = SOSCertificate(sp.Integer(2), (), (sp.Integer(1),), sp.ImmutableMatrix([[2]]))
    assert verify_sos_certificate(2, certificate)
