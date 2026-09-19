"""Bounded PSLQ candidate discovery with mandatory exact verification."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from math import gcd, lcm

import mpmath
import sympy as sp

from ._cost import ExactBudget, stage_allowed
from ._errors import EXACT_METHOD_ERRORS


@dataclass(frozen=True)
class PSLQRelation:
    """An integer relation that has passed an independent exact verifier."""

    terms: tuple[sp.Expr, ...]
    coeffs: tuple[int, ...]


def _primitive(values: list[int]) -> tuple[int, ...]:
    g = 0
    for value in values:
        g = gcd(g, abs(int(value)))
    if g == 0:
        return tuple(values)
    out = [int(v) // g for v in values]
    first = next((v for v in out if v), 1)
    if first < 0:
        out = [-v for v in out]
    return tuple(out)


def _linear_parts(term: sp.Expr, budget: ExactBudget):
    """Extract a small linear relation basis suitable for bounded PSLQ."""
    parts = term.args if term.is_Add else (term,)
    if len(parts) < 2 or len(parts) > budget.max_pslq_terms:
        return None
    bases = []
    coeffs = []
    for part in parts:
        coeff, base = part.as_coeff_Mul(rational=True)
        if not coeff.is_Rational:
            return None
        bases.append(base)
        coeffs.append(sp.Rational(coeff))
    scale = 1
    for coeff in coeffs:
        scale = lcm(scale, int(coeff.q))
    target = _primitive([int(coeff * scale) for coeff in coeffs])
    return tuple(bases), target


def pslq_zero_relation(term: sp.Expr,
                       verifier: Callable[[sp.Expr], bool | None],
                       budget: ExactBudget | None = None) -> PSLQRelation | None:
    """Find a candidate additive relation and require exact verification.

    PSLQ only selects a relation worth checking. The returned object exists
    only if ``verifier`` independently proves the corresponding exact symbolic
    relation to be zero.
    """
    budget = budget or ExactBudget()
    term = sp.sympify(term)
    if not stage_allowed(term, "pslq", budget):
        return None
    parsed = _linear_parts(term, budget)
    if parsed is None:
        return None
    bases, target = parsed
    old_dps = mpmath.mp.dps
    try:
        mpmath.mp.dps = budget.pslq_digits
        vals = [mpmath.mpf(str(sp.N(base, budget.pslq_digits))) for base in bases]
        rel = mpmath.pslq(mpmath.matrix(vals), tol=mpmath.mpf(10) ** (-(budget.pslq_digits - 15)),
                          maxcoeff=budget.pslq_coeff, maxsteps=200)
    except EXACT_METHOD_ERRORS:
        return None
    finally:
        mpmath.mp.dps = old_dps
    if rel is None:
        return None
    coeffs = _primitive([int(v) for v in rel])
    if coeffs != target:
        return None
    candidate = sp.Add(*(sp.Integer(c) * base for c, base in zip(coeffs, bases)), evaluate=False)
    try:
        if verifier(candidate) is True:
            return PSLQRelation(bases, coeffs)
    except EXACT_METHOD_ERRORS:
        return None
    return None
