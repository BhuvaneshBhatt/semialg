from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class GroebnerPrecondResult:
    input_polynomials: tuple[sp.Expr, ...]
    reduced_polynomials: tuple[sp.Expr, ...]
    changed: bool
    order: str
    variables: tuple[sp.Symbol, ...]


def groebner_precondition(
    polys: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    order: str = "lex",
) -> GroebnerPrecondResult:
    polys = tuple(sp.expand(p) for p in polys)
    vars_ = tuple(variables)
    if not polys or not vars_:
        return GroebnerPrecondResult(polys, polys, False, order, vars_)
    try:
        G = sp.groebner(polys, *vars_, order=order)
        reduced = tuple(sp.expand(g) for g in G.polys)
        return GroebnerPrecondResult(polys, reduced, reduced != polys, order, vars_)
    except Exception:
        return GroebnerPrecondResult(polys, polys, False, order, vars_)
