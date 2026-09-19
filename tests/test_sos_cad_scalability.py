import sympy as sp

from semialg.cad_algorithms.decomposition import decomp_collins_complete
from semialg.cad_algorithms.projection.collins import build_collins_proj_set
from semialg.formula import parse_formula
from semialg.partial.qe import lazy_find_inst_form
from semialg.sdp_backends import external_sdp_search, reconstruct_exact_gram
from semialg.sos_certificates import (
    SOSCertificate,
    plan_sos_search,
    sparse_sos_monomial_basis,
    verify_psd_exact,
    verify_sos_certificate,
)


def test_newton_basis_reduces_dense_gram_without_losing_certificate():
    x, y = sp.symbols("x y")
    polynomial = x**8 + y**8
    basis = sparse_sos_monomial_basis(polynomial, (x, y))
    plan = plan_sos_search(polynomial, (x, y))
    assert len(basis) == 5
    assert plan.gram_dimension == 5
    assert plan.dense_gram_dimension == 15
    assert plan.basis_method == "newton_polytope"
    gram = sp.diag(1, 0, 0, 0, 1)
    cert = SOSCertificate(polynomial, (x, y), basis, sp.ImmutableMatrix(gram))
    assert verify_sos_certificate(polynomial, cert, (x, y))


def test_exact_ldl_psd_handles_singular_and_indefinite_matrices():
    singular = sp.Matrix([[1, 1, 0], [1, 1, 0], [0, 0, 0]])
    indefinite = sp.Matrix([[0, 1], [1, 0]])
    result = verify_psd_exact(singular)
    assert result.verified
    assert result.method == "ldl_congruence"
    assert 0 in result.pivots
    assert not verify_psd_exact(indefinite).verified


def test_exact_gram_recovery_projects_numeric_hint_to_affine_identity():
    x, y = sp.symbols("x y")
    polynomial = x**4 + 2 * x**2 * y**2 + y**4
    basis = (x**2, x * y, y**2)
    numeric = [[1.00000001, 0.0, 0.99999999], [0.0, 0.0, 0.0], [0.99999999, 0.0, 1.00000001]]
    gram = reconstruct_exact_gram(polynomial, (x, y), basis, numeric, max_denominator=1000)
    assert gram is not None
    cert = SOSCertificate(polynomial, (x, y), basis, gram)
    assert verify_sos_certificate(polynomial, cert, (x, y))


def test_external_sdp_backend_declines_cleanly_without_cvxpy():
    x, y = sp.symbols("x y")
    result = external_sdp_search(x**4 + y**4, (x, y), backend="clarabel")
    assert not result.success
    assert result.status in {"cvxpy_unavailable", "solver_unavailable_or_failed"}


def test_collins_projection_uses_factors_instead_of_redundant_product():
    x, y = sp.symbols("x y")
    polynomial = (x - 1) * (x + 1) * (y - 1)
    tower = build_collins_proj_set((polynomial,), (x, y))
    top = {sp.factor(poly.as_expr()) for poly in tower.original_polynomials}
    assert top == {x - 1, x + 1, y - 1}


def test_lazy_existential_cad_stops_before_full_decomposition():
    x, y = sp.symbols("x y")
    formula = sp.And(x**2 < 1, y**2 < 1)
    lazy = lazy_find_inst_form((x, y), parse_formula(formula))
    full = decomp_collins_complete((x**2 - 1, y**2 - 1), lazy.stats.variables)
    assert lazy.found
    assert lazy.stats.stopped_early
    assert lazy.stats.evaluated_leaf_cells < len(full.cells_by_level[2])


def test_search_api_accepts_named_external_backend_and_preserves_exact_boundary():
    from semialg.sos_certificates import search_sos_certificate

    x, y = sp.symbols("x y")
    result = search_sos_certificate(x**4 + y**4, (x, y), backend="clarabel")
    assert not result.complete
    assert result.method == "sos_clarabel"
    assert any(note.startswith("backend_status:") for note in result.notes)
