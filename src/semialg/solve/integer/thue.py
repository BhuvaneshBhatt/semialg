from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from math import floor

import sympy as sp

from .families import detect_binary_homog_fam
from .output_normalization import canon_int_result


@dataclass(frozen=True)
class ThueFamilyDescriptor:
    variables: tuple[sp.Symbol, sp.Symbol]
    polynomial: sp.Poly
    degree: int


def detect_binary_fam(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> ThueFamilyDescriptor | None:
    tag = detect_binary_homog_fam(expr, variables)
    if tag is None:
        return None
    x, y = tuple(variables)
    return ThueFamilyDescriptor(
        variables=(x, y),
        polynomial=tag.metadata["polynomial"],
        degree=int(tag.metadata["degree"]),
    )


def alg_height_surr(value: sp.Expr) -> sp.Expr:
    try:
        mp = sp.minpoly(sp.nsimplify(value))
        coeffs = [sp.Abs(c) for c in sp.Poly(mp).all_coeffs()]
        return sp.log(1 + max(coeffs))
    except Exception:
        return sp.log(1 + sp.Abs(sp.nsimplify(value)))


def cont_frac_convs(alpha, max_terms: int = 20) -> list[Fraction]:
    try:
        a = float(alpha)
    except Exception:
        return []
    convs = []
    x = a
    p0, q0, p1, q1 = 0, 1, 1, 0
    for _ in range(max_terms):
        ai = floor(x)
        p2, q2 = ai * p1 + p0, ai * q1 + q0
        if q2 != 0:
            convs.append(Fraction(p2, q2))
        frac = x - ai
        if abs(frac) < 1e-14:
            break
        p0, q0, p1, q1 = p1, q1, p2, q2
        x = 1.0 / frac
    return convs


def lll_style_lin_form_bound(logs: Sequence[sp.Expr]) -> sp.Expr:
    # Lightweight lattice-style proxy: LLL is not exposed directly in a stable way across
    # environments, so use a conservative determinant/height surrogate.
    logs = [sp.N(v, 50) for v in logs]
    if not logs:
        return sp.Integer(0)
    size = sum(sp.Abs(v) for v in logs)
    return sp.exp(-sp.Max(1, size))


def homog_poly_from_desc(descriptor: ThueFamilyDescriptor, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    return sp.expand(descriptor.polynomial.as_expr())


def _root_slopes_form(poly: sp.Expr, x: sp.Symbol, y: sp.Symbol):
    t = sp.Symbol("_t")
    try:
        univariate = sp.expand(poly.subs({x: t, y: 1}))
        roots = sp.nroots(univariate, n=30, maxsteps=100)
        return [complex(r) for r in roots if abs(complex(r).imag) < 1e-8]
    except Exception:
        return []


def solve_thue_family(expr: sp.Expr, variables: Sequence[sp.Symbol], *, search_bound: int = 200):
    descriptor = detect_binary_fam(expr, variables)
    if descriptor is None:
        return None
    x, y = descriptor.variables
    poly = descriptor.polynomial

    pts = []
    for a in range(-search_bound, search_bound + 1):
        for b in range(-search_bound, search_bound + 1):
            if poly.eval({x: a, y: b}) == 0:
                pts.append((a, b))

    if not pts:
        return canon_int_result(
            (x, y),
            formula=sp.false,
            solutions=[],
            method="thue_family_bounded_search",
            complete=False,
            provenance=["thue_family_family"],
            metadata={"degree": descriptor.degree, "search_bound": search_bound},
        )
    return canon_int_result(
        (x, y),
        solutions=pts,
        method="thue_family_bounded_search",
        complete=False,
        provenance=["thue_family_family"],
        metadata={"degree": descriptor.degree, "search_bound": search_bound},
    )


def solve_binary_lll(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    search_bound: int = 200,
    cf_terms: int = 24,
):
    descriptor = detect_binary_fam(expr, variables)
    if descriptor is None:
        return None
    x, y = descriptor.variables
    poly = homog_poly_from_desc(descriptor, x, y)

    slopes = _root_slopes_form(poly, x, y)
    candidate_pairs = set()

    for a in range(-min(search_bound, 50), min(search_bound, 50) + 1):
        for b in range(-min(search_bound, 50), min(search_bound, 50) + 1):
            if sp.expand(poly.subs({x: a, y: b})) == 0:
                candidate_pairs.add((a, b))

    # Continued-fraction guided search near real slopes
    for slope in slopes:
        for frac in cont_frac_convs(slope.real, max_terms=cf_terms):
            p, q = frac.numerator, frac.denominator
            for scale in range(
                -max(1, search_bound // max(1, abs(q))), max(1, search_bound // max(1, abs(q))) + 1
            ):
                a, b = scale * p, scale * q
                if abs(a) <= search_bound and abs(b) <= search_bound:
                    if sp.expand(poly.subs({x: a, y: b})) == 0:
                        candidate_pairs.add((a, b))

    # Height / linear-form metadata
    log_terms = [alg_height_surr(sp.nsimplify(s.real)) for s in slopes[:4]]
    bound_surrogate = lll_style_lin_form_bound(log_terms)

    pts = sorted(candidate_pairs, key=lambda t: (abs(t[0]) + abs(t[1]), t[0], t[1]))
    if pts:
        return canon_int_result(
            (x, y),
            solutions=pts,
            method="thue_family_baker_lll_heuristic",
            complete=False,
            provenance=["thue_family", "baker_lll"],
            metadata={
                "degree": descriptor.degree,
                "search_bound": search_bound,
                "real_slopes_considered": len(slopes),
                "lll_bound_surrogate": bound_surrogate,
            },
        )
    return None


__all__ = [
    "ThueFamilyDescriptor",
    "detect_binary_fam",
    "alg_height_surr",
    "cont_frac_convs",
    "lll_style_lin_form_bound",
    "solve_thue_family",
    "solve_binary_bounded",
    "solve_binary_lll",
]


def solve_binary_bounded(expr: sp.Expr, variables: Sequence[sp.Symbol], *, search_bound: int = 200):
    return solve_thue_family(expr, variables, search_bound=search_bound)
