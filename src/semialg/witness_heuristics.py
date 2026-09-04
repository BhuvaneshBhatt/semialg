"""Cheap witness searches with exact post-verification.

These routines are deliberately one-sided: success is a certified witness;
failure is never interpreted as infeasibility.
"""

from __future__ import annotations

import random
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import product

import sympy as sp

from ._linear_relations import certified_sign
from .algebraic.roots import isolate_real_roots
from .algebraic.sample_points import choose_sector_sample


@dataclass(frozen=True)
class WitnessSearchResult:
    variables: tuple[sp.Symbol, ...]
    point: tuple[tuple[sp.Symbol, sp.Expr], ...] | None
    method: str
    certified: bool
    notes: tuple[str, ...] = ()

    @property
    def assignment(self) -> dict[sp.Symbol, sp.Expr] | None:
        return dict(self.point) if self.point is not None else None

    @property
    def found(self) -> bool:
        return self.point is not None and self.certified


def _certified_negative_point(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    values: Sequence[sp.Expr],
) -> tuple[tuple[sp.Symbol, sp.Expr], ...] | None:
    """Return a normalized assignment only when negativity is certified exactly."""
    vars_ = tuple(variables)
    point = tuple(values)
    if len(vars_) != len(point):
        raise ValueError("witness coordinate count must match variables")
    assignment = dict(zip(vars_, point, strict=True))
    if certified_sign(sp.expand(polynomial).subs(assignment)) != -1:
        return None
    return tuple((variable, sp.simplify(assignment[variable])) for variable in vars_)


def _univariate_negative_sample(poly: sp.Expr, variable: sp.Symbol) -> sp.Expr | None:
    """Return an exact point where a rational univariate polynomial is negative."""
    expr = sp.expand(poly)
    if certified_sign(expr.subs(variable, 0)) == -1:
        return sp.Integer(0)
    try:
        roots = tuple(isolate_real_roots(sp.Poly(expr, variable, domain=sp.QQ)))
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        roots = ()
    # Rational points between real roots suffice because the sign is constant
    # on every complementary interval.
    ordered = list(roots)
    rational_tests: list[sp.Expr] = []
    if ordered:
        try:
            rational_tests.append(choose_sector_sample(None, ordered[0]).value)
            for left_sample, right_sample in zip(ordered, ordered[1:], strict=False):
                rational_tests.append(choose_sector_sample(left_sample, right_sample).value)
            rational_tests.append(choose_sector_sample(ordered[-1], None).value)
        except (TypeError, ValueError):
            # Exact sector sampling can fail only when the root representation
            # cannot be refined enough.  A witness search is one-sided, so in
            # that case we simply leave this strategy unsuccessful.
            return None
    else:
        rational_tests.extend((sp.Integer(-1), sp.Integer(1)))
    for candidate in rational_tests:
        candidate = sp.Rational(candidate)
        if certified_sign(expr.subs(variable, candidate)) == -1:
            return candidate
    return None


def odd_degree_negative_witness(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
) -> WitnessSearchResult:
    """Find a certified negative point from an odd degree in one variable.

    If the degree in ``x_i`` is odd, choose exact integer values for the other
    variables until the leading coefficient is nonzero; the resulting
    univariate polynomial has opposite signs at the two infinities.  An exact
    univariate sign search then supplies a witness.
    """
    vars_ = tuple(variables)
    expr = sp.expand(polynomial)
    try:
        sp.Poly(expr, *vars_, domain=sp.QQ)
    except (sp.PolynomialError, ValueError, TypeError):
        return WitnessSearchResult(vars_, None, "odd_degree", False, ("non_rational_polynomial",))
    for variable in vars_:
        degree = sp.degree(expr, variable)
        if degree is None or int(degree) <= 0 or int(degree) % 2 == 0:
            continue
        leading = sp.Poly(expr, variable).LC()
        others = tuple(v for v in vars_ if v != variable)
        # Deterministic bounded lattice search.  Failure within this bounded
        # search is deliberately not interpreted as a mathematical conclusion.
        for radius in range(0, 8):
            candidates = [0] if radius == 0 else [-radius, radius]
            for values in product(candidates, repeat=len(others)):
                assignment = dict(zip(others, map(sp.Integer, values), strict=True))
                if sp.simplify(leading.subs(assignment)) == 0:
                    continue
                univariate = sp.expand(expr.subs(assignment))
                sample = _univariate_negative_sample(univariate, variable)
                if sample is None:
                    continue
                assignment[variable] = sample
                point = _certified_negative_point(expr, vars_, tuple(assignment[v] for v in vars_))
                if point is not None:
                    return WitnessSearchResult(vars_, point, "odd_degree", True)
    return WitnessSearchResult(vars_, None, "odd_degree", False)


def random_line_negative_witness(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    attempts: int = 8,
    seed: int = 1234,
    coefficient_bound: int = 16,
) -> WitnessSearchResult:
    """Search deterministic pseudo-random rational lines for a negative point.

    The search is heuristic, but every returned point is verified exactly.
    No conclusion is drawn when all attempted lines fail.
    """
    if attempts < 0:
        raise ValueError("attempts must be nonnegative")
    if coefficient_bound < 1:
        raise ValueError("coefficient_bound must be positive")
    vars_ = tuple(variables)
    expr = sp.expand(polynomial)
    try:
        sp.Poly(expr, *vars_, domain=sp.QQ)
    except (sp.PolynomialError, ValueError, TypeError):
        return WitnessSearchResult(vars_, None, "random_line", False, ("non_rational_polynomial",))
    rng = random.Random(seed)
    t = sp.Dummy("witness_line")
    for _ in range(max(0, attempts)):
        direction = tuple(rng.randint(-coefficient_bound, coefficient_bound) for _ in vars_)
        if not any(direction):
            continue
        line = {v: sp.Integer(a) * t for v, a in zip(vars_, direction, strict=True)}
        univariate = sp.expand(expr.subs(line))
        sample = _univariate_negative_sample(univariate, t)
        if sample is None:
            continue
        values = tuple(sp.Integer(a) * sample for a in direction)
        point = _certified_negative_point(expr, vars_, values)
        if point is not None:
            return WitnessSearchResult(vars_, point, "random_line", True)
    return WitnessSearchResult(vars_, None, "random_line", False)


def find_negative_witness_fast(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    random_lines: int = 8,
    seed: int = 1234,
) -> WitnessSearchResult:
    """Run cheap one-sided negative-witness strategies in deterministic order."""
    if random_lines < 0:
        raise ValueError("random_lines must be nonnegative")
    odd = odd_degree_negative_witness(polynomial, variables)
    if odd.found:
        return odd
    return random_line_negative_witness(polynomial, variables, attempts=random_lines, seed=seed)


__all__ = [
    "WitnessSearchResult",
    "find_negative_witness_fast",
    "odd_degree_negative_witness",
    "random_line_negative_witness",
]
