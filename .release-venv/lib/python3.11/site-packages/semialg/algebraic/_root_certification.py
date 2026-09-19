"""Construction and replay of exact polynomial root-interval certificates."""

from __future__ import annotations

import sympy as sp

from ._root_certificate_types import PolynomialRootIntervalCertificate
from .roots import (
    _count_roots_in_interval,
    _descartes_variations_interval,
    _exact_poly_eval,
    _flint_squarefree_part,
    _qq_univariate,
)


def certify_polynomial_root_interval(
    poly: sp.Poly | sp.Expr,
    left: sp.Rational | int,
    right: sp.Rational | int,
    *,
    var: sp.Symbol | None = None,
    include_left: bool = True,
    include_right: bool = True,
) -> PolynomialRootIntervalCertificate:
    """Certify the number of distinct real roots in a rational interval.

    For rational polynomials the cheap path uses exact Descartes variations on
    the square-free part.  Variation 0 or 1 is already decisive.  Sturm counting
    is reserved for ambiguous variation counts, matching the initial-isolation
    tree's certification discipline.  This routine accepts only
    polynomials; transcendental root solving is outside semialg.
    """

    univar = _qq_univariate(poly, var)
    a = sp.Rational(left)
    b = sp.Rational(right)
    if a > b:
        raise ValueError("left endpoint must not exceed right endpoint")
    squarefree = _flint_squarefree_part(univar)
    left_root = _exact_poly_eval(squarefree, a) == 0
    right_root = _exact_poly_eval(squarefree, b) == 0
    if a == b:
        count = int(left_root and include_left and include_right)
        return PolynomialRootIntervalCertificate(
            univar,
            a,
            b,
            count,
            left_root,
            right_root,
            "exact-point-evaluation",
            include_left,
            include_right,
            0,
        )

    variations = _descartes_variations_interval(squarefree, a, b)
    endpoint_count = int(left_root and include_left) + int(right_root and include_right)
    if variations <= 1:
        interior = variations
        method = "descartes"
    else:
        interior = _count_roots_in_interval(squarefree, a, b, exclude_left=True, exclude_right=True)
        method = "descartes+sturm"
    return PolynomialRootIntervalCertificate(
        univar,
        a,
        b,
        int(interior + endpoint_count),
        left_root,
        right_root,
        method,
        include_left,
        include_right,
        variations,
    )


def verify_polynomial_root_interval_certificate(
    certificate: PolynomialRootIntervalCertificate,
    *,
    polynomial: sp.Poly | sp.Expr | None = None,
    var: sp.Symbol | None = None,
    left: sp.Rational | int | None = None,
    right: sp.Rational | int | None = None,
    include_left: bool | None = None,
    include_right: bool | None = None,
) -> bool:
    """Replay an interval-root certificate and optionally bind it to expected input.

    The optional expected-input arguments let callers reject an otherwise valid
    certificate for a different polynomial, variable, interval, or endpoint
    convention. Omitting them preserves self-contained replay.
    """

    if not isinstance(certificate, PolynomialRootIntervalCertificate):
        return False
    try:
        if polynomial is not None:
            expected = _qq_univariate(polynomial, var)
            stored = _qq_univariate(certificate.polynomial)
            if expected.gens != stored.gens or expected.as_expr() != stored.as_expr():
                return False
        elif var is not None and certificate.polynomial.gens != (var,):
            return False
        if left is not None and certificate.left != sp.Rational(left):
            return False
        if right is not None and certificate.right != sp.Rational(right):
            return False
        if include_left is not None and certificate.include_left is not bool(include_left):
            return False
        if include_right is not None and certificate.include_right is not bool(include_right):
            return False
        replay = certify_polynomial_root_interval(
            certificate.polynomial,
            certificate.left,
            certificate.right,
            include_left=certificate.include_left,
            include_right=certificate.include_right,
        )
    except (ArithmeticError, TypeError, ValueError, sp.PolynomialError):
        return False
    return (
        replay.root_count == certificate.root_count
        and replay.left_is_root == certificate.left_is_root
        and replay.right_is_root == certificate.right_is_root
        and replay.method == certificate.method
        and replay.include_left == certificate.include_left
        and replay.include_right == certificate.include_right
        and replay.descartes_variations == certificate.descartes_variations
    )
