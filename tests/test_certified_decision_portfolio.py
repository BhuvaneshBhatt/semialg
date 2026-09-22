import pytest
import sympy as sp

from semialg import (
    find_instance,
    is_satisfiable,
    polynomial_nonnegative,
    reduce_formula,
    zeng_negative_point,
)
from semialg.formula import parse_quant_form_text
from semialg.sos_certificates import (
    SOSCertificate,
    plan_sos_search,
    search_sos_certificate,
    verify_sos_certificate,
)


def test_exact_sos_certificate_verifier_accepts_valid_gram_and_rejects_bad_psd():
    x, y = sp.symbols("x y", real=True)
    basis = (x, y)
    good = SOSCertificate(x**2 + y**2, (x, y), basis, sp.ImmutableMatrix.eye(2), "test")
    assert verify_sos_certificate(x**2 + y**2, good, (x, y))

    bad = SOSCertificate(x**2 - y**2, (x, y), basis, sp.ImmutableMatrix([[1, 0], [0, -1]]), "test")
    assert not verify_sos_certificate(x**2 - y**2, bad, (x, y))


def test_optional_search_backend_is_accepted_only_after_exact_verification():
    x, y = sp.symbols("x y", real=True)

    class ExternalCertificate:
        exact_residual_zero = True
        exact_psd_checks = (True,)
        gram_matrices = (sp.ImmutableMatrix.eye(2),)

    class Result:
        success = True
        monomial_basis = (x, y)
        certificate = ExternalCertificate()

    def searcher(polynomial, variables):
        return Result()

    searched = search_sos_certificate(x**2 + y**2, (x, y), backend=searcher)
    assert searched.certified
    result = polynomial_nonnegative(
        x**2 + y**2, (x, y), strategy="sos", sos_backend=searcher, return_result=True
    )
    assert result.decision is True
    assert result.backend == "sos"
    assert result.certificate is not None


def test_zeng_options_return_certified_decision():
    x = sp.Symbol("x", real=True)
    result = polynomial_nonnegative(
        x**2 + 1,
        (x,),
        strategy="zeng",
        return_result=True,
        random_lines=3,
        seed=77,
    )

    assert result.decision is True
    assert result.backend.startswith("zeng")


def test_single_public_nonnegative_api_validates_random_lines():
    x = sp.Symbol("x", real=True)
    with pytest.raises(ValueError, match="random_lines must be nonnegative"):
        polynomial_nonnegative(x**2 + 1, (x,), random_lines=-1)


def test_zeng_ars_nonnegative_case():
    x, y = sp.symbols("x y", real=True)
    polynomial = (x**2 + y**2 - 1) ** 2
    result = zeng_negative_point(polynomial, (x, y), random_lines=0)
    assert result.complete is True
    assert result.has_negative_point is False
    assert result.method == "zeng+ars"


def test_zeng_ars_negative_case():
    x, y = sp.symbols("x y", real=True)
    polynomial = (x**2 + y**2 - 1) ** 2 - 1
    result = zeng_negative_point(polynomial, (x, y), random_lines=0)
    assert result.complete is True
    assert result.has_negative_point is True
    assert result.method == "zeng+ars"
    assert result.assignment is not None
    assert sp.simplify(polynomial.subs(result.assignment)).is_negative is True


def test_is_satisfiable_uses_ars_for_positive_dimensional_pure_equalities():
    x, y = sp.symbols("x y", real=True)
    result = is_satisfiable(sp.Eq(x**2 + y**2, 1), (x, y), return_result=True)
    assert result.satisfiable is True
    assert "ars" in result.method
    assert result.witness is not None


def test_find_instance_uses_ars_for_positive_dimensional_pure_equalities():
    x, y = sp.symbols("x y", real=True)
    result = find_instance(sp.Eq(x**2 + y**2, 1), (x, y), return_result=True)
    assert result.first() is not None
    assert "ars" in result.method
    assert sp.simplify((x**2 + y**2 - 1).subs(result.first())) == 0


def test_existential_equality_uses_ars():
    parsed = parse_quant_form_text("exists x. exists y. x^2 + y^2 = 1")
    result = reduce_formula(parsed, strategy="auto", return_result=True)
    assert result.result is sp.true
    assert result.method.startswith("ars")
    assert "strategy_selection" not in result.metadata


def test_portfolio_falls_back_to_cad():
    x, y = sp.symbols("x y", real=True)
    # Noncoercive and not SOS; Zeng's current specialized fragment is incomplete.
    polynomial = x**2 * y**2 + x
    result = polynomial_nonnegative(polynomial, (x, y), sos_backend="none", return_result=True)
    assert result.decision is False
    assert result.backend in {"zeng+odd_degree", "zeng+random_line", "cad"}
    assert result.witness is not None


def test_sos_planner_skips_exactly_cheaper_univariate_and_quadratic_cases():
    x, y = sp.symbols("x y", real=True)
    univariate = plan_sos_search(x**8 + 1, (x,))
    assert univariate.applicable and not univariate.launch
    assert univariate.reason == "univariate_exact_preferred"

    quadratic = plan_sos_search(x**2 + y**2 + 1, (x, y))
    assert quadratic.applicable and not quadratic.launch
    assert quadratic.reason == "quadratic_exact_preferred"
    assert quadratic.gram_dimension == 3
    assert quadratic.gram_variables == 6


def test_sos_planner_uses_sparse_newton_basis_before_budgeting():
    x, y, z, w = sp.symbols("x y z w", real=True)
    plan = plan_sos_search((x**2 + y**2 + z**2 + w**2) ** 3, (x, y, z, w))
    assert plan.degree == 6
    assert plan.variable_count == 4
    assert plan.dense_gram_dimension == 35
    assert plan.gram_dimension == 20
    assert plan.gram_variables == 210
    assert plan.basis_method == "newton_polytope"
    assert plan.launch
    assert plan.reason == "within_auto_budget"


def test_auto_portfolio_respects_planner_route():
    x, y = sp.symbols("x y", real=True)
    calls = []

    def should_not_run(polynomial, variables):
        calls.append((polynomial, tuple(variables)))
        raise AssertionError("automatic SOS planner should have skipped the backend")

    result = polynomial_nonnegative(
        x**2 + y**2 + 1, (x, y), sos_backend=should_not_run, return_result=True
    )
    assert result.decision is True
    assert not calls
    first = result.attempts[0]
    assert first.backend == "sos"
    assert "planner:quadratic_exact_preferred" in first.notes
    assert "gram_dimension=3" in first.notes


def test_explicit_sos_strategy_bypasses_auto_cost_planner():
    x, y = sp.symbols("x y", real=True)
    calls = []

    class ExternalCertificate:
        exact_residual_zero = True
        exact_psd_checks = (True,)
        gram_matrices = (sp.ImmutableMatrix.eye(3),)

    class Result:
        success = True
        monomial_basis = (sp.Integer(1), x, y)
        certificate = ExternalCertificate()

    def search(polynomial, variables):
        calls.append((polynomial, tuple(variables)))
        return Result()

    result = polynomial_nonnegative(
        1 + x**2 + y**2, (x, y), strategy="sos", sos_backend=search, return_result=True
    )
    assert result.decision is True
    assert calls
    assert result.backend == "sos"


def test_odd_degree_is_marked_inapplicable_without_sdp_search():
    x, y = sp.symbols("x y", real=True)

    def should_not_run(*_args, **_kwargs):
        raise AssertionError("no SDP")

    searched = search_sos_certificate(x**3 + y**2, (x, y), backend=should_not_run, use_planner=True)
    assert not searched.certified
    assert searched.plan is not None
    assert searched.plan.applicable is False
    assert searched.plan.reason == "odd_total_degree"


def test_rich_external_certificate_is_preferred_over_flattened_gram():
    x, y = sp.symbols("x y", real=True)

    class ExternalCertificate:
        gram_matrices = [sp.ImmutableMatrix([[1]]), sp.ImmutableMatrix([[1]])]
        exact_residual_zero = True
        exact_psd_checks = [True, True]

    class Result:
        success = True
        monomial_basis = (x, y)
        certificate = ExternalCertificate()

    def searcher(polynomial, variables):
        return Result()

    searched = search_sos_certificate(x**2 + y**2, (x, y), backend=searcher)
    assert searched.certified
    assert searched.certificate is not None
    assert searched.certificate.source == "external:certificate"
    assert searched.certificate.gram_matrix == sp.ImmutableMatrix.eye(2)


def test_rich_external_certificate_metadata_is_not_trusted_as_proof():
    x, y = sp.symbols("x y", real=True)

    class ExternalCertificate:
        # Claims certification, but this Gram matrix proves x^2 - y^2 instead.
        gram_matrices = [sp.ImmutableMatrix([[1, 0], [0, -1]])]
        exact_residual_zero = True
        exact_psd_checks = [True]

    class Result:
        success = True
        monomial_basis = (x, y)
        certificate = ExternalCertificate()

    def searcher(polynomial, variables):
        return Result()

    searched = search_sos_certificate(x**2 + y**2, (x, y), backend=searcher)
    assert not searched.certified
    assert searched.certificate is None
    assert "candidate_not_exactly_verified" in searched.notes


def test_rich_external_certificate_with_unknown_psd_metadata_is_rejected():
    x, y = sp.symbols("x y", real=True)

    class ExternalCertificate:
        gram_matrices = [sp.ImmutableMatrix.eye(2)]
        exact_residual_zero = True
        exact_psd_checks = [None]

    class Result:
        success = True
        monomial_basis = (x, y)
        certificate = ExternalCertificate()

    def searcher(polynomial, variables):
        return Result()

    searched = search_sos_certificate(x**2 + y**2, (x, y), backend=searcher)
    assert not searched.certified
