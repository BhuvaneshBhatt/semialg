from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .structural_keys import symbol_identity_key


@dataclass(frozen=True)
class IncidenceComponent:
    variables: tuple[sp.Symbol, ...]
    polynomial_indices: tuple[int, ...]


@dataclass(frozen=True)
class IncidenceAnalysis:
    variables: tuple[sp.Symbol, ...]
    supports: tuple[frozenset[sp.Symbol], ...]
    components: tuple[IncidenceComponent, ...]
    degrees: tuple[tuple[sp.Symbol, int], ...]

    @property
    def decomposable(self) -> bool:
        return len(self.components) > 1


def analyze_incidence(
    polys: Sequence[sp.Expr | sp.Poly], variables: Sequence[sp.Symbol]
) -> IncidenceAnalysis:
    """Analyze polynomial-variable incidence, connected blocks, and degree burden."""
    vars_ = tuple(variables)
    supports: list[frozenset[sp.Symbol]] = []
    adjacency = {v: set() for v in vars_}
    degree_sum = {v: 0 for v in vars_}
    for raw in polys:
        expr = raw.as_expr() if isinstance(raw, sp.Poly) else sp.sympify(raw)
        support = frozenset(v for v in vars_ if v in expr.free_symbols)
        supports.append(support)
        for left in support:
            adjacency[left].update(s for s in support if s != left)
        if support:
            try:
                poly = sp.Poly(sp.expand(expr), *vars_)
            except (sp.PolynomialError, TypeError, ValueError):
                poly = None
            for var in support:
                try:
                    degree_sum[var] += int(poly.degree(var)) if poly is not None else 1
                except (TypeError, ValueError, sp.PolynomialError):
                    degree_sum[var] += 1

    components: list[IncidenceComponent] = []
    unseen = set(vars_)
    while unseen:
        seed = min(unseen, key=symbol_identity_key)
        stack = [seed]
        comp: set[sp.Symbol] = set()
        while stack:
            item = stack.pop()
            if item in comp:
                continue
            comp.add(item)
            stack.extend(adjacency[item] - comp)
        unseen.difference_update(comp)
        comp_vars = tuple(sorted(comp, key=symbol_identity_key))
        pidx = tuple(i for i, support in enumerate(supports) if support & comp)
        components.append(IncidenceComponent(comp_vars, pidx))
    components.sort(key=lambda c: tuple(symbol_identity_key(v) for v in c.variables))
    return IncidenceAnalysis(
        variables=vars_,
        supports=tuple(supports),
        components=tuple(components),
        degrees=tuple((v, degree_sum[v]) for v in vars_),
    )


def sparse_variable_order(
    polys: Sequence[sp.Expr | sp.Poly], variables: Sequence[sp.Symbol]
) -> tuple[sp.Symbol, ...]:
    """Return a deterministic low-incidence-width elimination order.

    This is a minimum-fill-style greedy heuristic over the variable incidence
    graph.  It never changes semantics; it only supplies a CAD ordering seed.
    """

    vars_ = tuple(variables)
    analysis = analyze_incidence(polys, vars_)
    neighbors = {v: set() for v in vars_}
    occurrence = {v: 0 for v in vars_}
    degree = dict(analysis.degrees)
    for support in analysis.supports:
        for v in support:
            occurrence[v] += 1
            neighbors[v].update(s for s in support if s != v)
    remaining = set(vars_)
    order: list[sp.Symbol] = []
    while remaining:

        def score(v: sp.Symbol):
            ns = neighbors[v] & remaining
            missing = sum(
                1
                for a in ns
                for b in ns
                if symbol_identity_key(a) < symbol_identity_key(b) and b not in neighbors[a]
            )
            return (missing, len(ns), occurrence[v], degree[v], symbol_identity_key(v))

        chosen = min(remaining, key=score)
        ns = list(neighbors[chosen] & remaining)
        for i, left in enumerate(ns):
            for right in ns[i + 1 :]:
                neighbors[left].add(right)
                neighbors[right].add(left)
        remaining.remove(chosen)
        order.append(chosen)
    return tuple(order)


def independent_polynomial_components(
    polys: Sequence[sp.Expr | sp.Poly], variables: Sequence[sp.Symbol]
) -> tuple[tuple[tuple[sp.Symbol, ...], tuple[sp.Expr, ...]], ...]:
    """Partition a polynomial family into independent variable components."""
    exprs = tuple(p.as_expr() if isinstance(p, sp.Poly) else sp.sympify(p) for p in polys)
    analysis = analyze_incidence(exprs, variables)
    out = []
    for component in analysis.components:
        part = tuple(exprs[i] for i in component.polynomial_indices)
        out.append((component.variables, part))
    return tuple(out)


__all__ = [
    "IncidenceAnalysis",
    "IncidenceComponent",
    "analyze_incidence",
    "independent_polynomial_components",
    "sparse_variable_order",
]


def decompose_conjunctive_formula(
    formula: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[tuple[sp.Expr, tuple[sp.Symbol, ...]], ...]:
    """Split a pure conjunction into independent variable-incidence blocks.

    Returns an empty tuple when no useful split exists. Parameter/constant atoms
    that involve none of ``variables`` are conservatively attached to the first
    block so their truth conditions are preserved.
    """

    vars_ = tuple(variables)
    atoms = tuple(formula.args) if isinstance(formula, sp.And) else (formula,)
    polys = []
    for atom in atoms:
        if not getattr(atom, "is_Relational", False):
            return tuple()
        polys.append(sp.expand(atom.lhs - atom.rhs))
    analysis = analyze_incidence(polys, vars_)
    active = [
        c for c in analysis.components if any(analysis.supports[i] for i in c.polynomial_indices)
    ]
    if len(active) <= 1:
        return tuple()
    assigned: set[int] = set()
    result: list[tuple[sp.Expr, tuple[sp.Symbol, ...]]] = []
    for component in active:
        idxs = [
            i
            for i, support in enumerate(analysis.supports)
            if support and support <= set(component.variables)
        ]
        if not idxs:
            continue
        assigned.update(idxs)
        result.append((sp.And(*(atoms[i] for i in idxs)), component.variables))
    leftovers = [atoms[i] for i in range(len(atoms)) if i not in assigned]
    if leftovers:
        if not result:
            return tuple()
        first_formula, first_vars = result[0]
        result[0] = (sp.And(first_formula, *leftovers), first_vars)
    return tuple(result) if len(result) > 1 else tuple()


__all__.append("decompose_conjunctive_formula")
