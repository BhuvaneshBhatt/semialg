"""Permanent semantic contracts for local/approach geometry."""

import sympy as sp

from semialg import (
    CorrelatedMapImageResult,
    LocalAlgebraicStrata,
    PointInClosureResult,
    ProofDiagnostics,
    angular_map_image,
    correlated_map_image,
    local_algebraic_strata,
    point_in_closure,
    structured_proof_diagnostics,
)


def test_closure_handles_boundary_interior_exterior_and_deleted_point():
    x = sp.Symbol("x", real=True)
    assert point_in_closure(sp.And(x > 0, x < 1), (0,), (x,))
    assert point_in_closure(sp.And(x > 0, x < 1), (sp.Rational(1, 2),), (x,))
    assert not point_in_closure(sp.And(x > 0, x < 1), (2,), (x,))
    assert point_in_closure(sp.And(x**2 < 1, sp.Ne(x, 0)), (0,), (x,))


def test_closure_is_invariant_under_equivalent_formula_and_coordinate_rename():
    x, z = sp.symbols("x z", real=True)
    first = sp.And(x > 0, x < 1)
    equivalent = sp.And(x * (1 - x) > 0, x > -2, x < 3)
    assert point_in_closure(first, (0,), (x,)) == point_in_closure(equivalent, (0,), (x,))
    assert point_in_closure(first, (0,), (x,)) == point_in_closure(
        first.xreplace({x: z}), (0,), (z,)
    )


def test_closure_commutes_with_invertible_affine_change_of_coordinates():
    x, z = sp.symbols("x z", real=True)
    source = sp.And(x > -1, x < 2, sp.Ne(x, 0))
    transformed = source.subs(x, (z - 3) / 2)
    for point in (-1, 0, 2, 3):
        expected = point_in_closure(source, (point,), (x,))
        assert point_in_closure(transformed, (2 * point + 3,), (z,)) == expected


def test_correlated_image_rejects_cartesian_range_overapproximation():
    t, u, v = sp.symbols("t u v", real=True)
    result = correlated_map_image((t, t**2), sp.And(t >= -1, t <= 1), (t,), image_variables=(u, v))
    assert sp.simplify(result.formula.subs({u: 0, v: 1})) is sp.false
    assert sp.simplify(result.formula.subs({u: -1, v: 1})) is sp.true
    assert result.diagnostics["correlated"] is True
    assert result.proof_trace[:2] == ("exact_graph", "existential_projection")


def test_correlated_image_is_equivariant_under_source_parameter_rename():
    t, s, u, v = sp.symbols("t s u v", real=True)
    a = correlated_map_image((t, t**2), sp.And(t >= -1, t <= 1), (t,), image_variables=(u, v))
    b = correlated_map_image((s, s**2), sp.And(s >= -1, s <= 1), (s,), image_variables=(u, v))
    samples = ((0, 0), (1, 1), (sp.Rational(1, 2), sp.Rational(1, 4)), (0, 1))
    for pu, pv in samples:
        assert sp.simplify(a.formula.subs({u: pu, v: pv})) == sp.simplify(
            b.formula.subs({u: pu, v: pv})
        )


def test_angular_image_sphere_constraint_changes_attainable_set():
    x, y, u, v = sp.symbols("x y u v", real=True)
    sphere = angular_map_image((x, y), (x, y), image_variables=(u, v), sphere=True)
    free = angular_map_image((x, y), (x, y), image_variables=(u, v), sphere=False)
    assert sp.simplify(sphere.formula.subs({u: 2, v: 0})) is sp.false
    assert sp.simplify(free.formula.subs({u: 2, v: 0})) is sp.true
    diagnostics = structured_proof_diagnostics(sphere)
    assert diagnostics.used_complete_qe
    assert "sphere_constraint" in diagnostics.steps


def test_local_strata_distinguish_regular_point_crossing_and_cusp():
    x, y = sp.symbols("x y", real=True)
    regular = local_algebraic_strata((y - x**2,), {x: 0, y: 0}, (x, y))
    crossing = local_algebraic_strata((x * y,), {x: 0, y: 0}, (x, y))
    cusp = local_algebraic_strata((y**2 - x**3,), {x: 0, y: 0}, (x, y))
    assert not regular.branch_geometry.singular
    assert crossing.branch_geometry.singular and len(crossing.branch_geometry.branches) == 2
    assert cusp.branch_geometry.singular
    assert all(
        sp.simplify(s.formula.subs({x: 0, y: 0})) is sp.true for s in cusp.incident_singular_strata
    )


def test_structured_diagnostics_is_stable_and_deduplicates_steps():
    class Result:
        method = "complete_qe_example"
        diagnostics = {"critical_values": (1,), "tag": "kept"}
        proof_trace = (
            "exact_graph",
            "existential_projection",
            "exact_graph",
            "critical_value_geometry",
        )

    d = structured_proof_diagnostics(Result())
    assert d.steps == ("exact_graph", "existential_projection", "critical_value_geometry")
    assert d.used_complete_qe and d.used_critical_values
    assert d.metadata["tag"] == "kept"


def test_limit_support_result_types_have_direct_construction_contracts():
    x, u = sp.symbols("x u", real=True)
    image = CorrelatedMapImageResult((x,), (x,), sp.true, (u,), sp.Eq(u, x))
    assert image.certified and image.method == "complete_qe_correlated_image"
    diagnostics = ProofDiagnostics("demo", ("demo",), False, False, {"stable": True})
    assert diagnostics.metadata["stable"] is True
    crossing = local_algebraic_strata((x * u,), {x: 0, u: 0}, (x, u))
    rebuilt = LocalAlgebraicStrata(
        crossing.branch_geometry,
        crossing.incident_singular_strata,
        crossing.incident_dimension_strata,
    )
    assert rebuilt.complete and rebuilt.branch_geometry == crossing.branch_geometry


def test_point_in_closure_result_direct_construction_contract():
    x = sp.Symbol("x", real=True)
    result = PointInClosureResult(True, x > 0, {x: 0}, (x,), x >= 0)
    assert bool(result) and result.closure == (x >= 0)
