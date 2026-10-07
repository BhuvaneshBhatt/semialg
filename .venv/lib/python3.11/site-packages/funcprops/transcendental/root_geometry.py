"""Certified root counting from derivative geometry on finite intervals."""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import sympy as sp

from .real_sign import robust_real_sign


@dataclass(frozen=True)
class MonotoneRootCount:
    expression: sp.Expr
    variable: sp.Symbol
    domain: sp.Interval
    count: int | None
    complete: bool
    critical_points: tuple = ()
    source: str = "derivative-sign-cells"


def _sign_at(expr, var, p, side):
    try:
        v = (
            sp.limit(expr, var, p, dir=side)
            if p in (sp.oo, -sp.oo)
            else sp.simplify(expr.subs(var, p))
        )
    except (TypeError, ValueError, NotImplementedError):
        return None
    return robust_real_sign(v)


def derivative_cell_root_count(expr, var, domain):
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if not isinstance(domain, sp.Interval):
        return None
    d = sp.diff(expr, var)
    try:
        roots = sp.solveset(sp.Eq(d, 0), var, domain=domain)
    except (TypeError, ValueError, NotImplementedError):
        return None
    if roots is sp.S.EmptySet:
        crit = ()
    elif isinstance(roots, sp.FiniteSet):
        crit = tuple(
            sorted(
                (r for r in roots if domain.contains(r) is sp.S.true),
                key=sp.default_sort_key,
            )
        )
    else:
        return None
    points = (domain.start, *crit, domain.end)
    count = 0
    # Each derivative-zero-free cell is monotone if derivative sign is certifiable.
    for a, b in itertools.pairwise(points):
        sample = (
            (a + b) / 2
            if a not in (-sp.oo, sp.oo) and b not in (-sp.oo, sp.oo)
            else (b - 1 if a is -sp.oo else a + 1)
        )
        ds = robust_real_sign(sp.simplify(d.subs(var, sample)))
        if ds not in (-1, 1):
            return None
        sa = _sign_at(expr, var, a, "+")
        sb = _sign_at(expr, var, b, "-")
        if sa is None or sb is None:
            return None
        if sa == 0:
            count += 1
        if sa * sb < 0:
            count += 1
    if points[-1] not in (-sp.oo, sp.oo) and _sign_at(expr, var, points[-1], "-") == 0:
        count += 1
    return MonotoneRootCount(expr, var, domain, count, True, crit)
