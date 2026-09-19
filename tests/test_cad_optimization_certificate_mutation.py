"""Domain-specific mutation checks for CAD and optimization proof objects."""

from __future__ import annotations

from dataclasses import replace

import sympy as sp

from semialg import replay_certificate, semialgebraic_minimize
from semialg.cad_algorithms.bounds import (
    CADBoundLevelVerification,
    CADCellBoundsCertificate,
    DelineabilityCertificate,
    RootOrderCertificate,
)


def test_cad_bounds_certificate_rejects_each_boolean_proof_mutation():
    x = sp.Symbol("x")
    level = CADBoundLevelVerification(1, x, True, True, True, True, True, True)
    cert = CADCellBoundsCertificate((1,), (level,), True, True)
    assert cert.verify()
    assert not replace(cert, certified=False).verify()
    assert not replace(cert, source_path_consistent=False).verify()
    for field in (
        "dependencies_valid",
        "section_verified",
        "root_order_verified",
        "repr_consistent",
        "sample_contained",
        "openness_valid",
    ):
        assert not replace(cert, levels=(replace(level, **{field: False}),)).verify()


def test_cad_root_certificates_reject_proof_mutations():
    x = sp.Symbol("x")
    delineability = DelineabilityCertificate(
        x**2 - 2, x, 0, sign_invariant=True, stack_order_verified=True, sample_root_verified=True
    )
    assert delineability.verify()
    for field in ("sign_invariant", "stack_order_verified", "sample_root_verified"):
        assert not replace(delineability, **{field: False}).verify()
    assert not replace(delineability, root_index=-1).verify()

    ordering = RootOrderCertificate(x, None, 0, 1, True, True)
    assert ordering.verify()
    assert not replace(ordering, adjacent=False).verify()
    assert not replace(ordering, order_verified=False).verify()


def test_optimization_result_replay_rejects_domain_specific_mutations():
    x = sp.Symbol("x")
    result = semialgebraic_minimize(x**2, x >= 1, (x,), return_result=True)
    assert replay_certificate(result).verified is True
    mutations = (
        replace(result, value=result.value + 1),
        replace(result, objective=result.objective + 1),
        replace(result, points=({x: result.points[0][x] + 1},)),
        replace(result, attained=not result.attained),
    )
    for mutated in mutations:
        assert replay_certificate(mutated).verified is False
