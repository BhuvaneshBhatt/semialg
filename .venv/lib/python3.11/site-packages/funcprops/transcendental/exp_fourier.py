from __future__ import annotations

import itertools
from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..exceptions import SYMBOLIC_QUERY_ERRORS
from .real_sign import robust_real_sign


@dataclass(frozen=True)
class ExponentialFourierCount:
    expression: sp.Basic
    variable: sp.Symbol
    domain: sp.Interval
    root_count: int | None
    complete: bool
    sequence: tuple[sp.Basic, ...] = ()
    critical_points: tuple[sp.Basic, ...] = ()
    evidence: tuple = ()
    source: str = "normalized-exponential-derivative-sequence"


def _affine_exp_term(term, var):
    """Return ``(rate, coefficient)`` for c*exp(rate*var+offset)."""
    term = sp.expand_power_exp(term)
    if not term.has(var):
        return S.Zero, term
    exps = tuple(term.atoms(sp.exp))
    if len(exps) != 1:
        return None
    e = exps[0]
    quotient = sp.simplify(term / e)
    if quotient.has(var):
        return None
    try:
        p = sp.Poly(sp.expand(e.args[0]), var)
    except SYMBOLIC_QUERY_ERRORS:
        return None
    if p.degree() > 1:
        return None
    rate = sp.simplify(p.nth(1))
    offset = sp.simplify(p.nth(0))
    coeff = sp.simplify(quotient * sp.exp(offset))
    if rate.has(var) or coeff.has(var):
        return None
    if rate.is_real is not True or coeff.is_real is False:
        return None
    return rate, coeff


def exponential_polynomial_terms(expr, var):
    """Parse a restricted real exponential polynomial.

    Supported inputs are finite sums of ``c*exp(a*x+b)`` and constants, with exact
    real coefficients/rates. Equal rates are combined.  This intentionally excludes
    polynomial prefactors and nested exponentials.
    """
    expr = sp.expand(sp.sympify(expr))
    terms = {}
    for term in sp.Add.make_args(expr):
        item = _affine_exp_term(term, var)
        if item is None:
            return None
        rate, coeff = item
        if coeff.is_Float:
            return None
        terms[rate] = sp.simplify(terms.get(rate, S.Zero) + coeff)
    terms = {a: c for a, c in terms.items() if robust_real_sign(c) != 0}
    if not terms:
        return ()
    rates = list(terms)
    try:
        rates.sort(key=lambda q: float(sp.N(q, 50)))
    except (TypeError, ValueError):
        return None
    for a, b in itertools.pairwise(rates):
        if robust_real_sign(b - a) != 1:
            return None
    return tuple((a, terms[a]) for a in rates)


def _normalized_derivative(expr, var):
    terms = exponential_polynomial_terms(expr, var)
    if terms is None or len(terms) <= 1:
        return None
    a0 = terms[0][0]
    normalized = sp.expand(sp.exp(-a0 * var) * expr)
    derivative = sp.simplify(sp.diff(normalized, var))
    return normalized, derivative


def _finite_exact_roots(expr, var, domain):
    try:
        roots = sp.solveset(sp.Eq(expr, 0), var, domain=domain)
    except SYMBOLIC_QUERY_ERRORS:
        return None
    if roots is S.EmptySet:
        return ()
    if not isinstance(roots, sp.FiniteSet):
        return None
    vals = []
    for r in roots:
        if r.is_real is False:
            return None
        contained = domain.contains(r)
        if contained is S.false:
            continue
        if contained is not S.true:
            return None
        vals.append(r)
    try:
        vals.sort(key=lambda q: float(sp.N(q, 50)))
    except (TypeError, ValueError):
        return None
    return tuple(vals)


def _point_sign(expr, var, point, direction=None):
    if point is S.NegativeInfinity or point is S.Infinity:
        try:
            value = sp.limit(expr, var, point, dir=direction or "+")
        except SYMBOLIC_QUERY_ERRORS:
            return None
    else:
        try:
            value = sp.simplify(expr.subs(var, point))
        except SYMBOLIC_QUERY_ERRORS:
            return None
    return robust_real_sign(value)


def exponential_fourier_count(expr, var, domain):
    """Certify the number of distinct real zeros of a restricted exp polynomial.

    The proof is recursive.  Multiplication by the positive factor exp(-a0*x)
    preserves zeros; differentiating the normalized expression removes its lowest
    exponential term.  Exact derivative roots partition the interval into cells on
    which the normalized expression is strictly monotone.  Endpoint/critical-point
    signs then count exactly one or zero roots per cell.  The method declines to
    certify whenever a derivative root or sign cannot be established exactly.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if not isinstance(domain, sp.Interval):
        return None
    terms = exponential_polynomial_terms(expr, var)
    if terms is None or not terms:
        return None
    if len(terms) == 1:
        sign = robust_real_sign(terms[0][1])
        if sign in (-1, 1):
            return ExponentialFourierCount(expr, var, domain, 0, True, (expr,))
        return None

    nd = _normalized_derivative(expr, var)
    if nd is None:
        return None
    normalized, derivative = nd
    derivative_roots = _finite_exact_roots(derivative, var, domain)
    if derivative_roots is None:
        return None
    derivative_roots = tuple(
        r for r in derivative_roots if r != domain.start and r != domain.end
    )

    # Recursively require the derivative sequence itself to be in the supported
    # class.  This is proof evidence that no derivative roots were missed.
    child = exponential_fourier_count(derivative, var, domain)
    if child is None or not child.complete or child.root_count != len(derivative_roots):
        return None

    points = [domain.start, *derivative_roots, domain.end]
    signs = []
    for point in points:
        if point == domain.start and domain.left_open:
            s = _point_sign(normalized, var, point, "+")
        elif point == domain.end and domain.right_open:
            s = _point_sign(normalized, var, point, "-")
        else:
            s = _point_sign(normalized, var, point)
        if s is None:
            return None
        signs.append(s)

    attained = [
        not domain.left_open,
        *([True] * len(derivative_roots)),
        not domain.right_open,
    ]
    count = sum(
        1
        for s, is_attained in zip(signs, attained, strict=True)
        if s == 0 and is_attained
    )
    # A critical-point root belongs to both neighboring closures; count it once and
    # do not also count a sign-change root through that endpoint.
    for left, right in itertools.pairwise(signs):
        if left and right and left != right:
            count += 1

    sequence = (normalized, *child.sequence)
    return ExponentialFourierCount(
        expr,
        var,
        domain,
        count,
        True,
        sequence,
        derivative_roots,
        evidence=(("signs", tuple(signs)), ("derivative", child)),
    )


__all__ = [
    "ExponentialFourierCount",
    "exponential_fourier_count",
    "exponential_polynomial_terms",
]
