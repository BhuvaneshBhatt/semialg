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


def test_zeng_unsupported_noncoercive_case_remains_incomplete_when_fast_search_disabled():
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
