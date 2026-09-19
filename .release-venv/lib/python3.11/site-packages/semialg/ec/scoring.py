from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ..structural_keys import symbol_identity_key


@dataclass(frozen=True)
class ECScore:
    expr: sp.Expr
    total_degree: int
    monomial_count: int
    symbol_count: int
    main_variable: sp.Symbol | None = None
    main_degree: int = 0
    projection_burden: int = 0

    @property
    def key(self) -> tuple[int, int, int, int, int, tuple[str, str], str]:
        main_key = (
            ("", "") if self.main_variable is None else symbol_identity_key(self.main_variable)
        )
        return (
            self.projection_burden,
            self.main_degree,
            self.total_degree,
            self.monomial_count,
            self.symbol_count,
            main_key,
            sp.srepr(self.expr),
        )


def score_eq_cons(expr: sp.Expr, variables: Sequence[sp.Symbol] = ()) -> ECScore:
    expanded = sp.expand(expr)
    explicit_gens = tuple(variables)
    gens = explicit_gens or tuple(sorted(expanded.free_symbols, key=symbol_identity_key))
    poly = sp.Poly(expanded, *gens) if gens else sp.Poly(expanded)
    main_variable = None
    main_degree = 0
    if gens:
        active = [(var, int(poly.degree(var))) for var in gens if poly.degree(var) > 0]
        if active:
            # CAD projects suffix variables first: prefer an EC whose highest
            # active variable has small degree and sparse coefficients.
            main_variable, main_degree = max(active, key=lambda item: gens.index(item[0]))
    coeff_count = 1
    if main_variable is not None:
        try:
            coeff_count = len(sp.Poly(expanded, main_variable).all_coeffs())
        except (sp.PolynomialError, TypeError, ValueError):
            coeff_count = len(poly.terms())
    projection_burden = (
        max(1, main_degree) * max(1, coeff_count) * max(1, len(poly.monoms()))
        if explicit_gens
        else 0
    )
    return ECScore(
        expr=expanded,
        total_degree=int(poly.total_degree()),
        monomial_count=len(poly.monoms()),
        symbol_count=len(expanded.free_symbols),
        main_variable=main_variable,
        main_degree=main_degree,
        projection_burden=projection_burden,
    )


def rank_eq_cons(
    exprs: Sequence[sp.Expr], variables: Sequence[sp.Symbol] = ()
) -> tuple[ECScore, ...]:
    uniq = tuple(dict.fromkeys(sp.expand(e) for e in exprs))
    return tuple(sorted((score_eq_cons(e, variables) for e in uniq), key=lambda s: s.key))


def choose_eq_cons(
    exprs: Sequence[sp.Expr], policy: str = "lowest_degree", variables: Sequence[sp.Symbol] = ()
) -> sp.Expr | None:
    ranked = rank_eq_cons(exprs, variables)
    if not ranked:
        return None
    if policy == "first":
        return sp.expand(exprs[0])
    if policy == "sparsest":
        return min(
            ranked,
            key=lambda s: (s.monomial_count, s.total_degree, s.symbol_count, sp.srepr(s.expr)),
        ).expr
    if policy in {"projection", "auto"}:
        return ranked[0].expr
    return min(
        ranked,
        key=lambda s: (s.total_degree, s.monomial_count, s.symbol_count, sp.srepr(s.expr)),
    ).expr


__all__ = ["ECScore", "score_eq_cons", "rank_eq_cons", "choose_eq_cons"]
