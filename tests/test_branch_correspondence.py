"""Ordered branch limits, complex births and exact Sturm germs."""

from dataclasses import replace

import pytest
import sympy as sp

from semialg.algebraic.samples import RationalSample
from semialg.topology.branch_correspondence import (
    certify_branch_correspondence,
    replay_branch_correspondence,
)


@pytest.mark.parametrize("source,index,expected", [(1, 0, 2), (-1, 0, -1), (-1, 1, -1), (-1, 2, 2)])
def test_complex_birth_does_not_redirect_branch(source, index, expected):
    x, y = sp.symbols("x y", real=True)
    certificate = certify_branch_correspondence(
        ((y + 1) ** 2 + x) * (y - 2),
        y,
        index,
        (x,),
        (RationalSample(sp.S(source)),),
        (sp.S.Zero,),
    )
    assert certificate.value == expected
    assert replay_branch_correspondence(certificate)
    assert not replay_branch_correspondence(
        replace(certificate, target_index=1 - certificate.target_index)
    )


@pytest.mark.parametrize("side,expected", [(1, -1), (-1, 1)])
def test_irreducible_collision_uses_sturm_germ(side, expected):
    x, y = sp.symbols("x y", real=True)
    for index in (0, 1):
        certificate = certify_branch_correspondence(
            (y**2 - 1) ** 2 + x * y,
            y,
            index,
            (x,),
            (RationalSample(sp.Rational(side, 4)),),
            (sp.S.Zero,),
        )
        assert certificate.value == expected
        assert certificate.continuation_ranks
        assert replay_branch_correspondence(certificate)
        assert not replay_branch_correspondence(replace(certificate, continuation_ranks=(99,)))


@pytest.mark.parametrize(
    "side,index,expected", [(1, 0, -sp.oo), (1, 1, 0), (-1, 0, 0), (-1, 1, sp.oo)]
)
def test_degree_loss(side, index, expected):
    x, y = sp.symbols("x y", real=True)
    proof = certify_branch_correspondence(
        x * y**2 + y, y, index, (x,), (RationalSample(sp.S(side)),), (sp.S.Zero,)
    )
    assert proof.value == expected
    assert replay_branch_correspondence(proof)
    assert not replay_branch_correspondence(replace(proof, target_index=2))


@pytest.mark.parametrize(
    "polynomial,count", [(lambda x, y: x * y * y - 1, 2), (lambda x, y: (x * y * y - 1) ** 2, 2)]
)
def test_both_infinite_ends(polynomial, count):
    x, y = sp.symbols("x y", real=True)
    proofs = [
        certify_branch_correspondence(
            polynomial(x, y), y, i, (x,), (RationalSample(sp.S.One),), (sp.S.Zero,)
        )
        for i in range(count)
    ]
    assert tuple(proof.value for proof in proofs) == (-sp.oo, sp.oo)
    assert all(replay_branch_correspondence(proof) for proof in proofs)


def test_persistent_multiple_source_branch():
    x, y = sp.symbols("x y", real=True)
    certificate = certify_branch_correspondence(
        (y**2 - x) ** 2, y, 1, (x,), (RationalSample(sp.S.One),), (sp.S.Zero,)
    )
    assert certificate.source_spectrum.multiplicities == (2, 2)
    assert certificate.value == 0


def test_collision_certificate_fields():
    from ._certificate_mutation import assert_certificate_fields_reject_mutation

    x, y = sp.symbols("x y", real=True)
    certificate = certify_branch_correspondence(
        (y**2 - 1) ** 2 + x * y, y, 0, (x,), (RationalSample(sp.Rational(1, 4)),), (sp.S.Zero,)
    )
    assert_certificate_fields_reject_mutation(certificate, replay_branch_correspondence)


@pytest.mark.parametrize(
    "kind", ["hyperbola", "two_branches", "upper_sector", "lower_sector", "quadratic_escape"]
)
def test_closure_at_infinity(kind):
    from semialg import cad

    x, y = sp.symbols("x y", real=True)
    fiber = {
        "hyperbola": sp.Eq(x * y, 1),
        "two_branches": sp.Eq(x * y * y + y, 0),
        "upper_sector": x * y > 1,
        "lower_sector": x * y < 1,
        "quadratic_escape": sp.Eq(x * y * y, 1),
    }[kind]
    result = cad(sp.And(x > 0, fiber), (x, y), operation="closure", return_result=True)
    function = result.as_function()
    for ordinate in (-2, -1, 0, 1, 2):
        expected = kind == "lower_sector" or (kind == "two_branches" and ordinate == 0)
        assert function(0, ordinate) is expected
        assert function(-1, ordinate) is False
    if kind == "upper_sector":
        assert function(1, 1) and not function(1, 0)
    if kind == "lower_sector":
        assert function(1, 1) and not function(1, 2)


@pytest.mark.parametrize("multiplicity", [1, 2])
@pytest.mark.parametrize("index,expected", [(0, -sp.oo), (1, sp.sqrt(3))])
def test_algebraic_degree_loss(multiplicity, index, expected):
    x, y = sp.symbols("x y", real=True)
    polynomial = ((x - sp.sqrt(2)) * y * y + y - sp.sqrt(3)) ** multiplicity
    proof = certify_branch_correspondence(
        polynomial,
        y,
        index,
        (x,),
        (RationalSample(sp.S(2)),),
        (sp.sqrt(2),),
    )
    assert proof.value == expected
    assert proof.source_spectrum.multiplicities == (multiplicity, multiplicity)
    assert replay_branch_correspondence(proof)


def test_unequal_persistent_multiplicities():
    x, y = sp.symbols("x y", real=True)
    polynomial = (x * y + 1) ** 2 * (y - 2) ** 3
    for index, expected in [(0, -sp.oo), (1, sp.S(2))]:
        proof = certify_branch_correspondence(
            polynomial,
            y,
            index,
            (x,),
            (RationalSample(sp.S.One),),
            (sp.S.Zero,),
        )
        assert proof.source_spectrum.multiplicities == (2, 3)
        assert proof.value == expected
        assert replay_branch_correspondence(proof)


@pytest.mark.parametrize("power", [1, 2])
def test_algebraic_parameter_content(power):
    x, y = sp.symbols("x y", real=True)
    proof = certify_branch_correspondence(
        (x - sp.sqrt(2)) ** power * (y - sp.sqrt(3)),
        y,
        0,
        (x,),
        (RationalSample(sp.S(2)),),
        (sp.sqrt(2),),
    )
    assert proof.value == sp.sqrt(3)
    assert sp.expand(proof.polynomial - y + sp.sqrt(3)) == 0
    assert replay_branch_correspondence(proof)
