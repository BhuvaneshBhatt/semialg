from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class ECScore:
    expr: sp.Expr
    total_degree: int
    monomial_count: int
    symbol_count: int

    @property
    def key(self) -> tuple[int, int, int, str]:
        return (self.total_degree, self.monomial_count, self.symbol_count, sp.srepr(self.expr))


def score_eq_cons(expr: sp.Expr) -> ECScore:
    poly = sp.Poly(sp.expand(expr))
    return ECScore(
        expr=sp.expand(expr),
        total_degree=poly.total_degree(),
        monomial_count=len(poly.monoms()),
        symbol_count=len(expr.free_symbols),
    )


def rank_eq_cons(exprs: Sequence[sp.Expr]) -> tuple[ECScore, ...]:
    uniq = tuple(dict.fromkeys(sp.expand(e) for e in exprs))
    return tuple(sorted((score_eq_cons(e) for e in uniq), key=lambda s: s.key))


def choose_eq_cons(exprs: Sequence[sp.Expr], policy: str = "lowest_degree") -> sp.Expr | None:
    ranked = rank_eq_cons(exprs)
    if not ranked:
        return None
    if policy == "first":
        return ranked[0].expr
    if policy == "sparsest":
        return min(
            ranked,
            key=lambda s: (s.monomial_count, s.total_degree, s.symbol_count, sp.srepr(s.expr)),
        ).expr
    return ranked[0].expr


__all__ = ["ECScore", "score_eq_cons", "rank_eq_cons", "choose_eq_cons"]
