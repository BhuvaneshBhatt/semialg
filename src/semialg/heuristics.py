from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence

import sympy as sp

from .formula import Formula, equational_constraints, formula_polynomials


def suggest_var_polys(
    polys: Sequence[sp.Expr],
    variables: Iterable[sp.Symbol] | None = None,
    *,
    strategy: str = "degree",
) -> tuple[sp.Symbol, ...]:
    polys = [sp.expand(poly) for poly in polys]
    if variables is None:
        vars_set = sorted(
            {sym for poly in polys for sym in poly.free_symbols}, key=lambda s: s.name
        )
    else:
        vars_set = list(variables)
    if strategy == "name":
        return tuple(sorted(vars_set, key=lambda s: s.name))

    degree_score = Counter()
    occurrence_score = Counter()
    for poly in polys:
        for sym in vars_set:
            if sym in poly.free_symbols:
                occurrence_score[sym] += 1
                try:
                    degree_score[sym] += sp.Poly(poly, *vars_set).degree(sym)
                except Exception:
                    degree_score[sym] += 1
    if strategy == "occurrence":
        ordered = sorted(
            vars_set, key=lambda sym: (occurrence_score[sym], degree_score[sym], sym.name)
        )
    else:
        ordered = sorted(
            vars_set, key=lambda sym: (degree_score[sym], occurrence_score[sym], sym.name)
        )
    return tuple(ordered)


def suggest_variable_order(formula: Formula, *, strategy: str = "degree") -> tuple[sp.Symbol, ...]:
    polys = formula_polynomials(formula)
    vars_ = set(sym for poly in polys for sym in poly.free_symbols)
    ecs = set(sym for expr in equational_constraints(formula) for sym in expr.free_symbols)
    base = list(suggest_var_polys(polys, sorted(vars_, key=lambda s: s.name), strategy=strategy))
    if strategy == "ec":
        base.sort(key=lambda sym: (sym not in ecs, base.index(sym)))
    return tuple(base)
