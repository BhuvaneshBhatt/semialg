from dataclasses import replace

import sympy as sp

from semialg.algebraic_decomposition import (
    EquidimensionalDecomposition,
    equidimensional_decomposition,
    verify_decomposition_certificate,
)
from semialg.region_analysis import _region_singular_locus_result, region_singular_locus


def _truth_at(formula, point):
    return sp.simplify(formula.subs(point))


def test_complete_decomposition_carries_replayable_certificate():
    x, y, z = sp.symbols("x y z", real=True)
    result = equidimensional_decomposition((x * y, x * z), (x, y, z))

    assert result.complete
    assert result.certificate is not None
    assert verify_decomposition_certificate(result.certificate)
    assert "monomial_minimal_primes" in result.certificate.methods


def test_incomplete_budget_limited_certificate_replays_with_same_budget():
    x, y, z = sp.symbols("x y z", real=True)
    result = equidimensional_decomposition((x * y, x * z), (x, y, z), max_pieces=1)

    assert not result.complete
    assert result.certificate is not None
    assert result.certificate.max_pieces == 1
    assert verify_decomposition_certificate(result.certificate)


def test_certificate_detects_tampering():
    x, y = sp.symbols("x y", real=True)
    result = equidimensional_decomposition((x * y,), (x, y))
    certificate = result.certificate
    assert certificate is not None

    tampered = replace(certificate, complete=not certificate.complete)
    assert not verify_decomposition_certificate(tampered)

    tampered_digest = replace(certificate, digest="0" * 64)
    assert not verify_decomposition_certificate(tampered_digest)


def test_region_analysis_refuses_complete_boolean_without_certificate():
    x, y = sp.symbols("x y", real=True)

    def uncertified_complete(_equations, variables, **_kwargs):
        return EquidimensionalDecomposition(
            variables=tuple(variables),
            pieces=(),
            complete=True,
            certificate=None,
        )

    result = _region_singular_locus_result(
        sp.Eq(x * y, 0), (x, y), decomposition_provider=uncertified_complete
    )
    assert not result.complete
    assert result.formula is None


def test_real_dimension_drop_union_has_only_true_component_intersection_singularity():
    x, y, z = sp.symbols("x y z", real=True)
    # Reduced real set: plane z=0 union line x=y=0.
    formula = sp.Eq((x**2 + y**2) * z, 0)
    singular = region_singular_locus(formula, (x, y, z))

    assert _truth_at(singular, {x: 0, y: 0, z: 0}) is sp.true
    assert _truth_at(singular, {x: 0, y: 0, z: 1}) is sp.false
    assert _truth_at(singular, {x: 1, y: 0, z: 0}) is sp.false


def test_nonreduced_real_dimension_drop_boundary_preserves_reduced_geometry():
    x, y, z = sp.symbols("x y z", real=True)
    formula = sp.And(sp.Eq((x**2 + y**2) ** 2 * z**3, 0), z >= 0)
    singular = region_singular_locus(formula, (x, y, z))

    assert _truth_at(singular, {x: 0, y: 0, z: 0}) is sp.true
    assert _truth_at(singular, {x: 0, y: 0, z: 2}) is sp.false
    assert _truth_at(singular, {x: 2, y: 0, z: 0}) is sp.false
