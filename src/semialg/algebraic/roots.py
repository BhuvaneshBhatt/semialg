from __future__ import annotations

from collections.abc import Iterable

import sympy as sp

from .cache import CACHE, poly_key
from .intervals import RationalInterval
from .samples import AlgebraicRoot


def _as_univariate_poly(poly: sp.Poly | sp.Expr, var: sp.Symbol | None = None) -> sp.Poly:
    if isinstance(poly, sp.Poly):
        if poly.is_multivariate:
            if var is None:
                raise ValueError("a variable is required for multivariate expressions")
            return sp.Poly(poly.as_expr(), var, domain="EX")
        return poly
    if var is None:
        symbols = sorted(poly.free_symbols, key=lambda s: s.name)
        if len(symbols) != 1:
            raise ValueError("a variable is required unless the expression is univariate")
        var = symbols[0]
    return sp.Poly(poly, var, domain="EX")


def _root_multiplicities(poly: sp.Poly) -> dict[sp.Expr, int]:
    out: dict[sp.Expr, int] = {}
    try:
        _, factors = sp.factor_list(poly.as_expr())
    except Exception:
        return out
    var = poly.gens[0]
    for factor_expr, mult in factors:
        factor_poly = sp.Poly(
            factor_expr, var, domain=poly.domain if poly.domain != sp.EX else "EX"
        )
        try:
            roots = sp.real_roots(factor_poly.as_expr())
        except Exception:
            roots = []
        for root in roots:
            out[root] = mult
    return out


def rational_intv_around(root: sp.Expr, *, digits: int = 80) -> RationalInterval:
    if root.is_Rational:
        return RationalInterval(sp.Rational(root), sp.Rational(root))
    value = sp.N(root, digits)
    center = sp.Rational(str(value))
    radius = sp.Rational(1, 10**20)
    return RationalInterval(center - radius, center + radius)


def _poly_intervals(poly: sp.Poly) -> list[tuple[RationalInterval, int]]:
    try:
        raw = poly.intervals(eps=sp.Rational(1, 10**30))
    except Exception:
        return []
    intervals: list[tuple[RationalInterval, int]] = []
    for bounds, mult in raw:
        left, right = bounds
        intervals.append((RationalInterval(sp.Rational(left), sp.Rational(right)), int(mult)))
    return intervals


def _dedupe_sorted_roots(roots: Iterable[sp.Expr]) -> list[sp.Expr]:
    ordered = sorted(roots, key=lambda root: sp.N(root, 100))
    unique: list[sp.Expr] = []
    for root in ordered:
        if not any(sp.simplify(root - old) == 0 for old in unique):
            unique.append(root)
    return unique


def isolate_real_roots(
    poly: sp.Poly | sp.Expr, var: sp.Symbol | None = None
) -> tuple[AlgebraicRoot, ...]:
    """Return explicit algebraic samples for the real roots of a univariate polynomial."""

    univar = _as_univariate_poly(poly, var)
    if univar.degree() <= 0:
        return ()
    key = (poly_key(univar), str(univar.gens[0]))
    CACHE.stats.calls += 1
    cached = CACHE.roots.get(key)
    if cached is not None:
        CACHE.stats.cache_hits += 1
        return cached  # type: ignore[return-value]
    CACHE.stats.cache_misses += 1

    expr = sp.expand(univar.as_expr())
    try:
        raw_roots = sp.real_roots(expr)
    except Exception:
        raw_roots = [
            sp.nsimplify(sp.re(root))
            for root in sp.nroots(univar)
            if abs(complex(sp.N(sp.im(root), 50))) < 1e-30
        ]
    roots = _dedupe_sorted_roots(raw_roots)
    intervals = _poly_intervals(univar)
    multiplicities = _root_multiplicities(univar)
    samples: list[AlgebraicRoot] = []
    for index, root in enumerate(roots):
        mult = next(
            (m for candidate, m in multiplicities.items() if sp.simplify(candidate - root) == 0), 1
        )
        interval = rational_intv_around(root)
        if index < len(intervals):
            interval, interval_mult = intervals[index]
            mult = interval_mult or mult
        samples.append(
            AlgebraicRoot(
                polynomial=univar,
                interval=interval,
                root_index=index,
                multiplicity=mult,
                root_expr=root,
            )
        )
    CACHE.roots[key] = tuple(samples)
    return tuple(samples)


def refine_isol_intv(root: AlgebraicRoot, *, steps: int = 4) -> AlgebraicRoot:
    """Return the same root with a narrower rational isolating interval.

    The refinement uses the exact algebraic expression to choose each half interval.
    """

    if root.interval.is_point():
        return root
    left, right = root.interval.left, root.interval.right
    expr = root.as_expr()
    for _ in range(max(steps, 0)):
        CACHE.stats.refinements += 1
        mid = sp.Rational(left + right, 2)
        if sp.N(expr - mid, 100) <= 0:
            right = mid
        else:
            left = mid
    return AlgebraicRoot(
        root.polynomial,
        RationalInterval(left, right),
        root.root_index,
        root.multiplicity,
        root.root_expr,
    )


def root_multiplicity(root: AlgebraicRoot) -> int:
    return root.multiplicity
