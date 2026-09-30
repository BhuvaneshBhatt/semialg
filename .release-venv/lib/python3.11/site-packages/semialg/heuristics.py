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


def _interaction_fill_score(
    polys: Sequence[sp.Expr], elimination_order: Sequence[sp.Symbol]
) -> tuple[int, int]:
    """Return fill edges and degree-weighted fill for a variable elimination order."""
    adjacency = {var: set() for var in elimination_order}
    for expr in polys:
        present = [var for var in elimination_order if var in expr.free_symbols]
        for index, left in enumerate(present):
            adjacency[left].update(present[index + 1 :])
            for right in present[index + 1 :]:
                adjacency[right].add(left)
    fill_edges = 0
    weighted_fill = 0
    remaining = set(elimination_order)
    for var in elimination_order:
        neighbors = sorted(adjacency[var] & remaining, key=symbol_identity_key)
        degree_weight = sum(_brown_key(polys, var)[:2])
        for i, left in enumerate(neighbors):
            for right in neighbors[i + 1 :]:
                if right not in adjacency[left]:
                    adjacency[left].add(right)
                    adjacency[right].add(left)
                    fill_edges += 1
                    weighted_fill += max(1, degree_weight)
        remaining.discard(var)
    return fill_edges, weighted_fill


def _abstract_projection_score(
    polys: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    lifting_order: Sequence[sp.Symbol],
) -> tuple[int, int]:
    """Estimate projection growth using multidegrees only."""
    variables = tuple(variables)
    positions = {var: index for index, var in enumerate(variables)}
    current = []
    for expr in polys:
        poly = sp.Poly(expr, *variables, domain="EX")
        current.append(tuple(max(0, int(poly.degree(var))) for var in variables))
    peak_count = len(current)
    total_degree_work = 0
    for variable in reversed(tuple(lifting_order)):
        position = positions[variable]
        active = [degrees for degrees in current if degrees[position] > 0]
        projected = [degrees for degrees in current if degrees[position] == 0]
        for degrees in active:
            eliminated_degree = degrees[position]
            projected.append(
                tuple(
                    0 if coordinate == position else degree * max(1, eliminated_degree)
                    for coordinate, degree in enumerate(degrees)
                )
            )
        for left_index, left in enumerate(active):
            for right in active[left_index + 1 :]:
                projected.append(
                    tuple(
                        0
                        if coordinate == position
                        else left[position] * right[coordinate] + right[position] * left[coordinate]
                        for coordinate in range(len(variables))
                    )
                )
        current = list(dict.fromkeys(projected))
        peak_count = max(peak_count, len(current))
        total_degree_work += sum(sum(degrees) for degrees in current)
    return peak_count, total_degree_work


def _cheap_auto_order(
    polys: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[tuple[sp.Symbol, ...], str]:
    """Choose between Brown and chordal using structural growth estimates."""
    brown = _brown_order(polys, variables)
    chordal = tuple(reversed(sparse_variable_order(polys, variables)))
    if chordal == brown or len(variables) < 3:
        return brown, "brown"

    brown_fill = _interaction_fill_score(polys, reversed(brown))
    chordal_fill = _interaction_fill_score(polys, reversed(chordal))
    if chordal_fill < brown_fill:
        return chordal, "chordal"
    if chordal_fill > brown_fill:
        return brown, "brown"

    brown_growth = _abstract_projection_score(polys, variables, brown)
    chordal_growth = _abstract_projection_score(polys, variables, chordal)
    return (chordal, "chordal") if chordal_growth < brown_growth else (brown, "brown")


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
        return _cheap_auto_order(polys, vars_set)[0]
    if strategy == "brown":
        return _brown_order(polys, vars_set)
    if strategy == "chordal":
        # sparse_variable_order is an elimination order; CAD APIs use lifting order.
        return tuple(reversed(sparse_variable_order(polys, vars_set)))

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
    """Suggest a CAD lifting order from structural or exact projection scores.

    ``strategy="auto"`` compares Brown and chordal candidates using interaction
    fill and an abstract multidegree projection-growth estimate.
    ``strategy="projection"`` scores every permutation up to
    ``exhaustive_limit`` and falls back to Brown above that limit.
    """

    polys = tuple(sp.expand(poly) for poly in polys)
    vars_ = tuple(variables)
    if strategy == "auto":
        order, selected = _cheap_auto_order(polys, vars_)
        score = tuple(sum(_brown_key(polys, v)[:3]) for v in reversed(order))
        return VariableOrderScore(order, selected, score)
    if strategy == "chordal":
        order = tuple(reversed(sparse_variable_order(polys, vars_)))
        score = tuple(sum(_brown_key(polys, v)[:3]) for v in reversed(order))
        return VariableOrderScore(order, "chordal", score)

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
    if strategy in {"auto", "brown", "chordal", "projection"}:
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
