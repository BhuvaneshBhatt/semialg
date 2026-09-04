from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from itertools import permutations

import sympy as sp

from .errors import InternalInvariantError
from .formula import Formula, equational_constraints, formula_polynomials
from .incidence import sparse_variable_order
from .structural_keys import symbol_identity_key


@dataclass(frozen=True)
class VariableOrderScore:
    order: tuple[sp.Symbol, ...]
    strategy: str
    score: tuple[int, ...]
    projection_polynomials: int | None = None
    sum_total_degree: int | None = None


def _brown_key(polys: Sequence[sp.Expr], var: sp.Symbol) -> tuple[int, int, int, str]:
    max_degree = 0
    max_term_degree = 0
    occurrences = 0
    for expr in polys:
        if var not in expr.free_symbols:
            continue
        occurrences += 1
        try:
            poly = sp.Poly(expr, *sorted(expr.free_symbols, key=symbol_identity_key), domain="EX")
            max_degree = max(max_degree, poly.degree(var))
            for monom, _coeff in poly.terms():
                if poly.gens:
                    idx = poly.gens.index(var) if var in poly.gens else None
                    if idx is not None and monom[idx]:
                        max_term_degree = max(max_term_degree, sum(monom))
        except (sp.PolynomialError, TypeError, ValueError):
            max_degree = max(max_degree, 1)
            max_term_degree = max(max_term_degree, 1)
    return (max_degree, max_term_degree, occurrences, var.name)


def _brown_order(polys: Sequence[sp.Expr], variables: Sequence[sp.Symbol]) -> tuple[sp.Symbol, ...]:
    # Collins eliminates the final variable first. Brown's rule chooses the
    # cheapest elimination variable first, hence reverse the ascending scores
    # into lifting order.
    eliminated = sorted(tuple(variables), key=lambda v: _brown_key(polys, v))
    return tuple(reversed(eliminated))


def _projection_score(polys: Sequence[sp.Expr], order: tuple[sp.Symbol, ...]) -> tuple[int, int]:
    from .cad_algorithms.projection.collins import build_collins_proj_set

    tower = build_collins_proj_set(polys, order)
    total_count = sum(len(level.polynomials) for level in tower.levels)
    sotd = 0
    for level in tower.levels:
        for poly in level.polynomials:
            try:
                sotd += sum(sum(monom) for monom, _ in poly.terms())
            except (TypeError, ValueError):
                sotd += max(0, poly.total_degree())
    return total_count, sotd


def rank_variables_by_polynomials(
    polys: Sequence[sp.Expr],
    variables: Iterable[sp.Symbol] | None = None,
    *,
    strategy: str = "degree",
) -> tuple[sp.Symbol, ...]:
    """Order variables using inexpensive polynomial incidence and degree data.

    ``strategy`` selects name, sparse-incidence, Brown, occurrence, or degree
    ordering. The routine does not build a projection set; use
    :func:`suggest_cad_variable_order` when projection cost should be measured.
    """
    polys = [sp.expand(poly) for poly in polys]
    if variables is None:
        vars_set = sorted(
            {sym for poly in polys for sym in poly.free_symbols}, key=symbol_identity_key
        )
    else:
        vars_set = list(variables)
    if strategy == "name":
        return tuple(sorted(vars_set, key=symbol_identity_key))
    if strategy == "auto":
        sparse = sparse_variable_order(polys, vars_set)
        brown = _brown_order(polys, vars_set)

        # Prefer the sparse seed when it begins with a lower-incidence variable;
        # otherwise retain Brown's established degree-first behavior.
        def occurrence(v):
            return sum(v in p.free_symbols for p in polys)

        return (
            sparse
            if (occurrence(sparse[0]), len(sparse)) < (occurrence(brown[0]), len(brown))
            else brown
        )
    if strategy == "brown":
        return _brown_order(polys, vars_set)

    degree_score = Counter()
    occurrence_score = Counter()
    for poly in polys:
        for sym in vars_set:
            if sym in poly.free_symbols:
                occurrence_score[sym] += 1
                try:
                    degree_score[sym] += sp.Poly(poly, *vars_set).degree(sym)
                except (sp.PolynomialError, TypeError, ValueError):
                    degree_score[sym] += 1
    if strategy == "occurrence":
        ordered = sorted(
            vars_set,
            key=lambda sym: (occurrence_score[sym], degree_score[sym], symbol_identity_key(sym)),
        )
    else:
        ordered = sorted(
            vars_set,
            key=lambda sym: (degree_score[sym], occurrence_score[sym], symbol_identity_key(sym)),
        )
    return tuple(ordered)


def suggest_cad_variable_order(
    polys: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    strategy: str = "auto",
    exhaustive_limit: int = 4,
) -> VariableOrderScore:
    """Suggest a CAD lifting order using Brown or exact projection-set scoring.

    ``strategy='auto'`` uses Brown's inexpensive heuristic. ``strategy='projection'``
    scores every permutation up to ``exhaustive_limit`` and falls back to Brown
    above that limit.
    """

    polys = tuple(sp.expand(poly) for poly in polys)
    vars_ = tuple(variables)
    use_projection = strategy == "projection"
    if use_projection and len(vars_) <= exhaustive_limit and len(vars_) > 1:
        best_order = vars_
        best_score: tuple[int, int] | None = None
        for order in permutations(vars_):
            score = _projection_score(polys, tuple(order))
            if best_score is None or score < best_score:
                best_order = tuple(order)
                best_score = score
        if best_score is None:
            raise InternalInvariantError("projection-order scoring produced no candidate")
        return VariableOrderScore(
            best_order,
            "projection",
            best_score,
            projection_polynomials=best_score[0],
            sum_total_degree=best_score[1],
        )
    order = _brown_order(polys, vars_)
    brown_score = tuple(sum(_brown_key(polys, v)[:3]) for v in reversed(order))
    return VariableOrderScore(order, "brown", brown_score)


def suggest_variable_order(formula: Formula, *, strategy: str = "degree") -> tuple[sp.Symbol, ...]:
    """Suggest a CAD variable order using the requested structural heuristic."""
    polys = formula_polynomials(formula)
    vars_ = set(sym for poly in polys for sym in poly.free_symbols)
    ecs = set(sym for expr in equational_constraints(formula) for sym in expr.free_symbols)
    if strategy in {"auto", "brown", "projection"}:
        return suggest_cad_variable_order(
            polys, sorted(vars_, key=symbol_identity_key), strategy=strategy
        ).order
    base = list(
        rank_variables_by_polynomials(
            polys, sorted(vars_, key=symbol_identity_key), strategy=strategy
        )
    )
    if strategy == "ec":
        base.sort(key=lambda sym: (sym not in ecs, base.index(sym)))
    return tuple(base)


__all__ = [
    "VariableOrderScore",
    "suggest_cad_variable_order",
    "rank_variables_by_polynomials",
    "suggest_variable_order",
]
