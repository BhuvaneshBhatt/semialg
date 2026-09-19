import sympy as sp

from semialg import find_negative_point, polynomial_nonnegative, zeng_negative_point
from semialg.reasoning import prove_nonnegative


def test_zeng_coercive_critical_values_certify_nonnegative():
    x, y = sp.symbols("x y")
    result = zeng_negative_point(x**2 + y**2 + 1, (x, y), random_lines=0)
    assert result.complete is True
    assert result.has_negative_point is False
    assert polynomial_nonnegative(x**2 + y**2 + 1, (x, y)) is True


def test_zeng_returns_exact_negative_critical_or_fast_witness():
    x, y = sp.symbols("x y")
    point = find_negative_point(x**2 + y**2 - 1, (x, y))
    assert point is not None
    assert sp.simplify((x**2 + y**2 - 1).subs(point)).is_negative is True


def test_zeng_odd_degree_fast_path_is_one_sided_and_exact():
    x, y = sp.symbols("x y")
    result = zeng_negative_point(x**3 + y**2, (x, y))
    assert result.complete is True
    assert result.has_negative_point is True
    assert result.assignment is not None
    assert sp.simplify((x**3 + y**2).subs(result.assignment)) < 0


def test_zeng_declines_noncoercive_case():
    x, y = sp.symbols("x y")
    result = zeng_negative_point(x**2 * y**2, (x, y), random_lines=0)
    assert result.has_negative_point is None
    assert result.complete is False


def test_global_sign_reasoning_can_consume_zeng_backend():
    x, y = sp.symbols("x y")
    # Not a syntactic sum of squares due to the cross term.
    assert prove_nonnegative(x**4 + y**4 + x * y + 3, (x, y)) is True


def test_coercivity_certificate_rejects_mixed_odd_exponent_leading_terms():
    x, y = sp.symbols("x y")
    polynomial = x**4 + y**4 + 10 * x**3 * y
    result = zeng_negative_point(polynomial, (x, y), random_lines=0)
    assert result.complete is False
    assert result.has_negative_point is None
    assert result.nonnegative is None
    assert sp.expand(polynomial.subs(y, -x)) == -8 * x**4


def test_coercivity_sphere_certificate_accepts_positive_mixed_leading_form():
    from semialg.polynomial_positivity import _coercivity_certificate

    x, y = sp.symbols("x y", real=True)
    leading = x**4 - x**2 * y**2 + y**4
    certificate = _coercivity_certificate(leading + x + 7, (x, y))
    assert certificate.certified is True
    assert certificate.method == "leading_form_sphere"
    assert certificate.leading_form == leading
    assert any(note.startswith("sphere_method=") for note in certificate.notes)


def test_zeng_uses_stronger_sphere_coercivity_certificate():
    x, y = sp.symbols("x y", real=True)
    polynomial = x**4 - x**2 * y**2 + y**4 + 1
    result = zeng_negative_point(polynomial, (x, y), random_lines=0)
    assert result.complete is True
    assert result.has_negative_point is False
    assert "coercivity_method=leading_form_sphere" in result.notes


def test_coercivity_sos_margin_route_precedes_sphere():
    from types import SimpleNamespace

    from semialg.polynomial_positivity import _coercivity_certificate

    x, y = sp.symbols("x y", real=True)
    leading = sp.expand((x**2 + x * y + y**2) ** 2)
    calls = []

    def fake_search(polynomial, variables, **kwargs):
        calls.append((sp.expand(polynomial), tuple(variables), kwargs))
        return SimpleNamespace(certified=True)

    certificate = _coercivity_certificate(leading + 1, (x, y), sos_searcher=fake_search)
    assert certificate.certified is True
    assert certificate.method == "leading_form_sos_margin"
    assert certificate.margin is not None and certificate.margin > 0
    assert calls


def test_coercivity_direction_disproof_avoids_sphere_for_degenerate_form():
    from semialg.polynomial_positivity import _coercivity_certificate

    x, y = sp.symbols("x y", real=True)
    certificate = _coercivity_certificate(x**2 * y**2 + 1, (x, y))
    assert certificate.certified is False
    assert certificate.method == "leading_form_direction"
    assert any(note.startswith("nonpositive_direction=") for note in certificate.notes)
