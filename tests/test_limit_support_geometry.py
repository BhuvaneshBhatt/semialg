import sympy as sp

from semialg import (
    CorrelatedMapImageResult,
    LocalAlgebraicStrata,
    ProofDiagnostics,
    angular_map_image,
    correlated_map_image,
    local_algebraic_strata,
    structured_proof_diagnostics,
)


def test_correlated_map_image_preserves_circle_relation():
    t = sp.symbols("t", real=True)
    u, v = sp.symbols("u v", real=True)
    r = correlated_map_image((t, t**2), sp.And(t >= -1, t <= 1), (t,), image_variables=(u, v))
    assert r.certified
    assert sp.simplify(r.formula.subs({u: sp.Rational(1, 2), v: sp.Rational(1, 4)})) is sp.true
    assert sp.simplify(r.formula.subs({u: sp.Rational(1, 2), v: sp.Rational(1, 3)})) is sp.false


def test_angular_map_image_keeps_coordinate_correlation():
    x, y = sp.symbols("x y", real=True)
    u, v = sp.symbols("u v", real=True)
    r = angular_map_image((x, y), (x, y), image_variables=(u, v))
    assert r.certified
    assert sp.simplify(r.formula.subs({u: 1, v: 0})) is sp.true
    assert sp.simplify(r.formula.subs({u: 1, v: 1})) is sp.false
    d = structured_proof_diagnostics(r)
    assert "sphere_constraint" in d.steps
    assert d.used_complete_qe


def test_local_algebraic_strata_records_crossing_branches():
    x, y = sp.symbols("x y", real=True)
    r = local_algebraic_strata((x * y,), {x: 0, y: 0}, (x, y))
    assert r.complete and r.branch_geometry.singular
    assert len(r.branch_geometry.branches) == 2
    assert r.branch_geometry.intersection is not None


def test_limit_support_structured_result_construction():
    x, u = sp.symbols("x u", real=True)
    image = CorrelatedMapImageResult((x,), (x,), sp.true, (u,), sp.Eq(u, x))
    assert image.formula == sp.Eq(u, x)
    proof = ProofDiagnostics("manual", ("manual",), False, False, {})
    assert proof.steps == ("manual",)
    strata = local_algebraic_strata((x * u,), {x: 0, u: 0}, (x, u))
    copied = LocalAlgebraicStrata(
        strata.branch_geometry, strata.incident_singular_strata, strata.incident_dimension_strata
    )
    assert copied.complete
