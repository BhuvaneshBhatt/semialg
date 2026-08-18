from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from ..validation.solution_checking import form_sat_by_assign


@dataclass(frozen=True)
class UnivarDecompWit:
    variable: sp.Symbol
    witness: object | None
    method: str


def _cand_values_comps(expr: sp.Expr, variable: sp.Symbol, domain) -> list[object]:
    candidates = []
    try:
        solset = sp.solveset(sp.Eq(expr, 0), variable, domain=domain)
        if isinstance(solset, sp.FiniteSet):
            candidates.extend(list(solset))
    except Exception:
        pass
    try:
        poly = sp.Poly(sp.expand(expr), variable)
        if poly.total_degree() > 0:
            candidates.extend(poly.all_roots())
    except Exception:
        pass
    out = []
    for c in candidates:
        if c not in out:
            out.append(c)
    return out


def find_univar_decomp_wit(
    equation_expr: sp.Expr, side_condition: sp.Expr, variable: sp.Symbol, *, domain=sp.Reals
) -> UnivarDecompWit:
    for candidate in _cand_values_comps(equation_expr, variable, domain):
        if form_sat_by_assign(
            sp.And(sp.Eq(equation_expr, 0), side_condition),
            {variable: candidate},
            domain=domain,
            check_numeric_equalities=True,
        ):
            return UnivarDecompWit(variable, candidate, "decomposition_search")
    return UnivarDecompWit(variable, None, "decomposition_search")


__all__ = ["UnivarDecompWit", "find_univar_decomp_wit"]
