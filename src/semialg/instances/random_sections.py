from __future__ import annotations

import random
from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class RandomSectionWitness:
    assignment: dict[sp.Symbol, object] | None
    attempts: int
    seed: int


def find_random_section_wit(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    domain=sp.Reals,
    seed: int = 7,
    attempts: int = 7,
) -> RandomSectionWitness:
    variables = tuple(variables)
    if len(variables) <= 1:
        return RandomSectionWitness(None, 0, seed)
    rng = random.Random(seed)
    section_symbol = sp.Symbol("_section_parameter", real=(domain == sp.Reals))
    scale = 8
    for attempt in range(1, attempts + 1):
        replacement = {var: rng.randint(-scale, scale) * section_symbol for var in variables}
        reduced = formula.subs(replacement)
        pieces = []
        if isinstance(reduced, sp.And):
            for atom in reduced.args:
                if isinstance(atom, sp.Equality):
                    pieces.append(atom.lhs - atom.rhs)
                else:
                    pieces.append(atom)
        elif isinstance(reduced, sp.Equality):
            pieces = [reduced.lhs - reduced.rhs]
        else:
            pieces = [reduced]
        try:
            sols = sp.solve(pieces, (section_symbol,), dict=True)
        except Exception:
            sols = []
        if isinstance(sols, list) and sols:
            value = sols[0].get(section_symbol)
            if value is not None:
                assignment = {
                    var: sp.simplify(expr.subs(section_symbol, value))
                    for var, expr in replacement.items()
                }
                return RandomSectionWitness(assignment, attempt, seed)
        scale *= 2
    return RandomSectionWitness(None, attempts, seed)


__all__ = ["RandomSectionWitness", "find_random_section_wit"]
