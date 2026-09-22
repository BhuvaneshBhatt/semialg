import sympy as sp

import semialg
from semialg.constraint_certificates import (
    component_constraint_descriptions,
    implied_polynomial_inequality,
    nonnegative_combination_certificate,
    redundant_polynomial_inequalities,
    verify_nonnegative_combination_certificate,
)
from semialg.sos_certificates import SOSCertificate


def test_nonnegative_combination():
    x, y = sp.symbols("x y", real=True)
    c = nonnegative_combination_certificate(2 * x + 3 * y, (x, y), (x, y), allow_sos=False)
    assert c is not None and verify_nonnegative_combination_certificate(c, (x, y))


def test_implied_and_redundant():
    x = sp.symbols("x", real=True)
    p = implied_polynomial_inequality(sp.And(x >= 0, x <= 1), 2 * x >= 0, (x,))
    assert p is not None and p.certified
    r = redundant_polynomial_inequalities(sp.And(x >= 0, 2 * x >= 0, x <= 1), (x,))
    assert any(item.constraint == (2 * x >= 0) for item in r)


def test_root_exports():
    assert semialg.implied_polynomial_inequality is implied_polynomial_inequality


def test_supplied_sos_certificate_and_components():
    x, y = sp.symbols("x y", real=True)
    sos = SOSCertificate(x**2 + y**2, (x, y), (x, y), sp.ImmutableMatrix.eye(2))
    c = nonnegative_combination_certificate(x**2 + y**2, (), (x, y), sos_certificate=sos)
    assert c is not None and c.sos_certificate == sos
    parts = component_constraint_descriptions(sp.And(sp.Eq(x * y, 0), x >= 0), (x, y))
    assert len(parts) == 2


def test_certificate_verifiers_reject_wrong_payload_types():
    x = sp.symbols("x", real=True)
    assert verify_nonnegative_combination_certificate(object(), (x,)) is False
    from semialg.sos_certificates import verify_sos_certificate

    assert verify_sos_certificate(x**2, object(), (x,)) is False


def test_sos_verifier_rejects_malformed_typed_payloads():
    from dataclasses import replace

    from semialg.sos_certificates import verify_sos_certificate

    x = sp.symbols("x", real=True)
    cert = SOSCertificate(x**2, (x,), (x,), sp.ImmutableMatrix([[1]]))
    assert verify_sos_certificate(x**2, replace(cert, variables=None), (x,)) is False
    assert verify_sos_certificate(x**2, replace(cert, monomial_basis=None), (x,)) is False
    assert verify_sos_certificate(x**2, replace(cert, gram_matrix="bad"), (x,)) is False
    assert verify_sos_certificate(x**2, replace(cert, polynomial=object()), (x,)) is False


def test_nonnegative_combination_verifier_rejects_malformed_typed_payloads():
    from dataclasses import replace

    from semialg.constraint_certificates import NonnegativeCombinationCertificate

    x = sp.symbols("x", real=True)
    cert = NonnegativeCombinationCertificate(x, (x,), (sp.Integer(1),), sp.Integer(0))
    assert verify_nonnegative_combination_certificate(replace(cert, premises=None), (x,)) is False
    assert (
        verify_nonnegative_combination_certificate(replace(cert, coefficients=None), (x,)) is False
    )
    assert verify_nonnegative_combination_certificate(replace(cert, target=object()), (x,)) is False
    assert (
        verify_nonnegative_combination_certificate(replace(cert, residual=object()), (x,)) is False
    )
