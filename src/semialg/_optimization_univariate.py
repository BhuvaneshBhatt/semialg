"""Complete polynomial calculus on exact univariate semialgebraic domains."""

from functools import cmp_to_key

import sympy as sp

from .exact_arithmetic import compare_exact_reals, exact_sign
from .interval_decomposition import (
    finite_real_roots,
    relational_polynomials,
    sample_between,
    truth_at,
)
from .optimization_results import OptimizationResult


def univariate_polynomial_optimum(objective, condition, variable, *, kind):
    """Certify extrema by exhaustive root isolation and exact polynomial limits.

    Every constraint root is a domain cut; its truth value retains isolated
    feasible points. Between cuts all constraint signs are constant. All
    derivative roots and finite/infinite endpoint limits exhaust possible
    extrema on each feasible interval. Any unavailable root or comparison
    causes this specialization to decline, never to claim partial coverage.
    """
    try:
        expressions = (objective, *relational_polynomials(condition))
    except TypeError:
        return None
    polynomials = []
    for expression in expressions:
        try:
            poly = sp.Poly(expression, variable, extension=True)
        except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed):
            return None
        if any(
            c.free_symbols or c.is_real is not True or c.is_algebraic is not True
            for c in poly.coeffs()
        ):
            return None
        polynomials.append(poly)
    polynomial = polynomials[0]
    try:
        cuts = tuple(
            sorted(
                set(
                    root
                    for poly in polynomials[1:]
                    for root in finite_real_roots(poly.as_expr(), variable)
                ),
                key=cmp_to_key(compare_exact_reals),
            )
        )
        ordered_cuts = []
        for cut in cuts:
            if not ordered_cuts or compare_exact_reals(ordered_cuts[-1], cut) != 0:
                ordered_cuts.append(cut)
        cuts = tuple(ordered_cuts)
        derivative = polynomial.diff()
        stationary = finite_real_roots(derivative.as_expr(), variable)
        candidates = []
        from .optimization import _best_candidates, _Candidate

        for point in cuts:
            if truth_at(condition, {variable: point}):
                candidates.append(
                    _Candidate(sp.simplify(polynomial.eval(point)), {variable: point}, True)
                )
        degree = polynomial.degree()

        def endpoint_value(point):
            if point not in (sp.oo, -sp.oo):
                return sp.simplify(polynomial.eval(point))
            if degree <= 0:
                return polynomial.as_expr()
            sign = exact_sign(polynomial.LC())
            if point == -sp.oo and degree % 2:
                sign = -sign
            return sign * sp.oo

        boundaries = (-sp.oo, *cuts, sp.oo)
        for left, right in zip(boundaries, boundaries[1:], strict=False):
            sample = sample_between(left, right)
            if not truth_at(condition, {variable: sample}):
                continue
            candidates.extend(
                (
                    _Candidate(endpoint_value(left), None, False),
                    _Candidate(endpoint_value(right), None, False),
                )
            )
            for root in stationary:
                if (left == -sp.oo or compare_exact_reals(left, root) < 0) and (
                    right == sp.oo or compare_exact_reals(root, right) < 0
                ):
                    candidates.append(
                        _Candidate(sp.simplify(polynomial.eval(root)), {variable: root}, True)
                    )
            if derivative.is_zero:
                candidates.append(_Candidate(polynomial.as_expr(), {variable: sample}, True))
        if candidates:
            value, points, attained = _best_candidates(candidates, kind=kind)
    except (NotImplementedError, sp.PolynomialError, TypeError, ValueError, ArithmeticError):
        return None
    if not candidates:
        raise ValueError("optimization domain is empty")
    return OptimizationResult(
        objective,
        (variable,),
        value,
        points,
        attained,
        kind,
        "exact_univariate_polynomial_calculus",
        {
            "constraint_formula": condition,
            "global_certificate": "complete_real_root_coverage_and_polynomial_limits",
            "domain_cut_count": len(cuts),
            "stationary_root_count": len(stationary),
        },
        True,
    )
