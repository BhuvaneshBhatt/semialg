"""Deterministic expression-sensitive proof routing and bounded complexity budgets."""

from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp


@dataclass
class ProofBudget:
    max_cost: int = 100
    spent: int = 0
    attempted: list[str] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)

    def allow(self, route, cost):
        if cost < 0:
            raise ValueError("cost must be nonnegative")
        if self.spent + cost > self.max_cost:
            self.skipped.append(route)
            return False
        self.spent += cost
        self.attempted.append(route)
        return True


def expression_cost(expr):
    e = sp.sympify(expr)
    ops = int(sp.count_ops(e, visual=False)) + 1
    len(e.free_symbols)
    len(e.args) if isinstance(e, sp.Piecewise) else 0
    try:
        syms = tuple(e.free_symbols)
        if syms and e.is_polynomial(*syms):
            max(0, int(sp.Poly(e, *syms).total_degree()))
    except (sp.PolynomialError, TypeError, ValueError):
        pass
    sum(1 for a in sp.preorder_traversal(e) if getattr(a, "is_Function", False))
    return ops


def route_proofs(expr, routes, *, budget=None):
    budget = budget or ProofBudget()
    accepted = []
    base = expression_cost(expr)
    for name, multiplier in routes:
        if budget.allow(name, max(1, int(base * multiplier))):
            accepted.append(name)
    return tuple(accepted), budget
