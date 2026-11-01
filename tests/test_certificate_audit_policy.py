"""Coverage-policy failures remain test failures, not rejected-proof successes."""

from dataclasses import dataclass

import pytest
import sympy as sp

from semialg.sos_certificates import SOSCertificate, verify_sos_certificate

from ._certificate_mutation import _policy, assert_certificate_fields_reject_mutation


def test_unregistered_fields():
    @dataclass
    class Unregistered:
        new_proof_field: int = 1

    with pytest.raises(AssertionError, match="unclassified"):
        _policy(Unregistered())


def test_unexpected_replay_error():
    x = sp.Symbol("x")
    certificate = SOSCertificate(x * x + 1, (x,), (sp.Integer(1), x), sp.ImmutableMatrix.eye(2))

    def verifier(candidate):
        if candidate != certificate:
            raise RuntimeError("unexpected verifier bug")
        return verify_sos_certificate(x * x + 1, candidate)

    with pytest.raises(RuntimeError, match="unexpected verifier bug"):
        assert_certificate_fields_reject_mutation(certificate, verifier)
