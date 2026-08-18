from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from itertools import permutations

import sympy as sp

from ..formula import Formula, formula_polynomials
from .features import ProblemFeatures


@dataclass(frozen=True)
class OrderScore:
    order: tuple[sp.Symbol, ...]
    score: int
    reason: str = ""


def _poly(poly: sp.Expr, order: Sequence[sp.Symbol]) -> sp.Poly | None:
    try:
        return sp.Poly(sp.expand(poly), *order)
    except Exception:
        return None


def _degree_in(poly: sp.Expr, sym: sp.Symbol) -> int:
    try:
        return int(sp.Poly(sp.expand(poly), sym).degree())
    except Exception:
        return 0


def _total_degree(poly: sp.Expr, order: Sequence[sp.Symbol]) -> int:
    pobj = _poly(poly, order)
    return int(pobj.total_degree()) if pobj is not None else 0


def brown_variable_order(
    polys: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Symbol, ...]:
    """Brown-style CAD variable ordering heuristic.

    Brown's heuristic eliminates variables using low maximum degree, then low
    total degree in terms containing that variable, then low occurrence count.
    The returned order is the corresponding CAD variable order.
    """

    vars_tuple = tuple(variables)

    def key(sym: sp.Symbol):
        max_degree = max((_degree_in(poly, sym) for poly in polys), default=0)
        degree_sum = sum(_degree_in(poly, sym) for poly in polys)
        occurrence = sum(1 for poly in polys if sym in poly.free_symbols)
        return (max_degree, degree_sum, occurrence, sym.name)

    return tuple(sorted(vars_tuple, key=key))


def sotd_score(order: Sequence[sp.Symbol], polys: Sequence[sp.Expr]) -> int:
    """Approximate sum-of-total-degrees score for an order.

    This is intentionally cheap: it scores the input family under the proposed
    order without constructing a full projection tower. It is suitable for
    ranking candidate orders before CAD construction.
    """

    return sum(_total_degree(poly, order) for poly in polys)


def ndrr_score(order: Sequence[sp.Symbol], polys: Sequence[sp.Expr]) -> int:
    """Cheap NDRR proxy: count distinct univariate real roots after projection to first variable.

    The full NDRR heuristic requires projection. This proxy keeps decomposition selection
    lightweight by looking at input polynomials that are univariate in the first
    variable under the proposed order.
    """

    if not order:
        return 0
    first = order[0]
    count = 0
    seen: set[str] = set()
    for poly in polys:
        if poly.free_symbols and poly.free_symbols <= {first}:
            try:
                roots = sp.Poly(poly, first).real_roots()
            except Exception:
                roots = ()
            for root in roots:
                key = sp.sstr(root)
                if key not in seen:
                    seen.add(key)
                    count += 1
    return count


def score_variable_order(
    order: Sequence[sp.Symbol], polys: Sequence[sp.Expr], *, ec_vars: Iterable[sp.Symbol] = ()
) -> OrderScore:
    order = tuple(order)
    ec_vars = set(ec_vars)
    position = {sym: i for i, sym in enumerate(order)}
    score = 0
    # Lower score is better. The components approximate Brown, SOTD, NDRR, and
    # EC placement without constructing a complete projection for every order.
    brown_rank = {sym: idx for idx, sym in enumerate(brown_variable_order(polys, order))}
    score += 50 * sum(abs(position[sym] - brown_rank[sym]) for sym in order)
    score += 10 * sotd_score(order, polys)
    score += 25 * ndrr_score(order, polys)
    for poly in polys:
        syms = [sym for sym in order if sym in poly.free_symbols]
        if not syms:
            continue
        width = position[syms[-1]] - position[syms[0]]
        score += width
        score += len(syms)
    # EC variables earlier in the CAD order tend to preserve lower-dimensional
    # structure for formula output and reduced projection selection.
    score += sum(position[sym] for sym in ec_vars if sym in position)
    return OrderScore(
        order=order, score=score, reason="Brown/SOTD/NDRR/EC composite; smaller is better"
    )


def cand_var_orders(
    features: ProblemFeatures,
    polys: Sequence[sp.Expr],
    *,
    limit: int = 12,
) -> tuple[OrderScore, ...]:
    vars_ = tuple(features.variables)
    if len(vars_) <= 1:
        return (OrderScore(order=vars_, score=0, reason="single variable"),)
    ec_vars = {sym for poly in polys for sym in poly.free_symbols} if features.has_ecs else set()
    sorted_vars = tuple(sorted(vars_, key=lambda s: s.name))
    candidates: set[tuple[sp.Symbol, ...]] = {
        tuple(vars_),
        sorted_vars,
        tuple(reversed(sorted_vars)),
        brown_variable_order(polys, vars_),
        tuple(reversed(brown_variable_order(polys, vars_))),
        tuple(
            sorted(
                vars_, key=lambda sym: (sum(sym in poly.free_symbols for poly in polys), sym.name)
            )
        ),
        tuple(
            sorted(
                vars_, key=lambda sym: (-sum(sym in poly.free_symbols for poly in polys), sym.name)
            )
        ),
    }
    if len(vars_) <= 5:
        all_orders = list(permutations(vars_))
        scored = sorted(
            (score_variable_order(order, polys, ec_vars=ec_vars) for order in all_orders),
            key=lambda item: item.score,
        )
        for item in scored[:limit]:
            candidates.add(item.order)
    scores = [score_variable_order(cand, polys, ec_vars=ec_vars) for cand in candidates]
    scores.sort(key=lambda s: (s.score, tuple(sym.name for sym in s.order)))
    return tuple(scores[:limit])


def choose_best_var_order(
    features: ProblemFeatures, polys: Sequence[sp.Expr]
) -> tuple[sp.Symbol, ...]:
    return cand_var_orders(features, polys)[0].order


def choose_best_form(formula: Formula, features: ProblemFeatures) -> tuple[sp.Symbol, ...]:
    return choose_best_var_order(features, tuple(formula_polynomials(formula)))


__all__ = [
    "OrderScore",
    "brown_variable_order",
    "cand_var_orders",
    "choose_best_form",
    "choose_best_var_order",
    "ndrr_score",
    "score_variable_order",
    "sotd_score",
]
