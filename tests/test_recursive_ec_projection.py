"""Recursive EC projection with degree and order conditions."""

from dataclasses import replace

import pytest
import sympy as sp

from semialg.cad_algorithms.projection.ec_hypotheses import replay_ec_projection
from semialg.cad_algorithms.reduced import decomp_form_reduced_safe
from semialg.formula import parse_formula
from semialg.qe.complete import qe_from_cad


@pytest.mark.parametrize("backend", ["mccallum", "lazard"])
def test_intermediate_non_ec_discriminant(backend):
    u, x, y, z = sp.symbols("u x y z", real=True)
    matrix = parse_formula(sp.And(sp.Eq(z - y, 0), sp.Eq(y - x, 0), y**2 > u))
    result = decomp_form_reduced_safe(matrix, (u, x, y, z), backend=backend)
    assert result.certificate.checked_refinement_cells == 0
    plan = next(plan for plan in result.cad.tower.metadata["ec_levels"] if plan.level == 3)
    assert plan.operator == "ec-order"
    assert u in plan.order_polynomials
    projected = qe_from_cad(result.cad, (u, x, y, z), (("exists", y), ("exists", z)), matrix)
    for uv in (-1, 0, 1, 4):
        for xv in (-2, 0, 1, 2):
            assert bool(projected.subs({u: uv, x: xv})) is (xv * xv > uv)
    altered = replace(plan, degree_coefficients=())
    plans = tuple(
        altered if item.level == 3 else item for item in result.cad.tower.metadata["ec_levels"]
    )
    tower = replace(
        result.cad.tower, metadata={**dict(result.cad.tower.metadata), "ec_levels": plans}
    )
    assert not replay_ec_projection(tower, matrix)


def test_resultant_implication_is_replayed():
    x, y, z = sp.symbols("x y z", real=True)
    matrix = parse_formula(
        sp.And(sp.Eq(x + y**2 + z, 0), sp.Eq(x - y**2 + z, 0), x**2 + y**2 + z**2 >= 1)
    )
    result = decomp_form_reduced_safe(matrix, (x, y, z), backend="mccallum")
    assert result.certificate.checked_refinement_cells == 0
    proofs = result.cad.tower.metadata["ec_implications"]
    derived = next(proof for proof in proofs if proof.polynomial == y)
    assert derived.eliminated_variable == z
    assert replay_ec_projection(result.cad.tower, matrix)
    projected = qe_from_cad(result.cad, (x, y, z), (("exists", y), ("exists", z)), matrix)
    for value in (-2, -1, sp.Rational(-1, 2), 0, sp.Rational(1, 2), 1, 2):
        assert bool(projected.subs(x, value)) is bool(2 * value**2 >= 1)
    altered = tuple(
        replace(proof, polynomial=y - 1) if proof == derived else proof for proof in proofs
    )
    tower = replace(
        result.cad.tower, metadata={**dict(result.cad.tower.metadata), "ec_implications": altered}
    )
    assert not replay_ec_projection(tower, matrix)


def test_resultant_product_does_not_imply_each_factor():
    x, y, z = sp.symbols("x y z", real=True)
    matrix = parse_formula(sp.And(sp.Eq(y * z, 0), x >= 0))
    result = decomp_form_reduced_safe(matrix, (x, y, z), backend="mccallum")
    selected = [plan.equational_constraint for plan in result.cad.tower.metadata["ec_levels"]]
    assert y not in selected
    assert z not in selected
    projected = qe_from_cad(result.cad, (x, y, z), (("exists", y), ("exists", z)), matrix)
    for value in (-1, 0, 1):
        assert bool(projected.subs(x, value)) is (value >= 0)


def test_hypothesis_field_mutations():
    from ._certificate_mutation import assert_certificate_fields_reject_mutation

    u, x, y, z = sp.symbols("u x y z", real=True)
    matrix = parse_formula(sp.And(sp.Eq(z - y, 0), sp.Eq(y - x, 0), y * y > u))
    source = decomp_form_reduced_safe(matrix, (u, x, y, z), backend="mccallum").cad.tower
    plans = source.metadata["ec_levels"]
    selected = next(i for i, plan in enumerate(plans) if plan.operator == "ec-order")

    def replay_plan(candidate):
        changed = (*plans[:selected], candidate, *plans[selected + 1 :])
        return replay_ec_projection(
            replace(source, metadata={**dict(source.metadata), "ec_levels": changed}), matrix
        )

    assert_certificate_fields_reject_mutation(plans[selected], replay_plan)
    implications = source.metadata["ec_implications"]

    def replay_implication(candidate):
        changed = (candidate, *implications[1:])
        return replay_ec_projection(
            replace(source, metadata={**dict(source.metadata), "ec_implications": changed}), matrix
        )

    assert_certificate_fields_reject_mutation(implications[0], replay_implication)


def test_resultant_implication_mutations():
    from ._certificate_mutation import assert_certificate_fields_reject_mutation

    x, y, z = sp.symbols("x y z", real=True)
    matrix = parse_formula(sp.And(sp.Eq(z * z - y, 0), sp.Eq(z * z + y, 0), x * x >= 1))
    source = decomp_form_reduced_safe(matrix, (x, y, z), backend="mccallum").cad.tower
    implications = source.metadata["ec_implications"]
    selected = next(
        i for i, claim in enumerate(implications) if claim.eliminated_variable is not None
    )

    def replay(candidate):
        changed = (*implications[:selected], candidate, *implications[selected + 1 :])
        return replay_ec_projection(
            replace(source, metadata={**dict(source.metadata), "ec_implications": changed}), matrix
        )

    assert_certificate_fields_reject_mutation(implications[selected], replay)
