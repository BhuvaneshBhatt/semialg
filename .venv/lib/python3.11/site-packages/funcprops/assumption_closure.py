"""Cheap, certified closure of assumptions used by property rules."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from .assumption_normalization import normalize_assumptions


@dataclass(frozen=True)
class AssumptionClosure:
    original: sp.Expr
    facts: frozenset[sp.Expr]
    condition: sp.Expr

    def entails(self, predicate):
        from .assumptions import entails

        p = sp.sympify(predicate)
        if p in self.facts:
            return True
        return entails(p, self.condition)


def _consequences(atom):
    out = {atom}
    if isinstance(atom, (sp.StrictGreaterThan, sp.StrictLessThan)):
        x, y = atom.lhs, atom.rhs
        if y == 0:
            out.add(sp.Ne(x, 0))
            out.add(sp.Contains(x, S.Reals))
            out.add(
                sp.Ge(x, 0) if isinstance(atom, sp.StrictGreaterThan) else sp.Le(x, 0)
            )
        elif x == 0:
            out.add(sp.Ne(y, 0))
            out.add(sp.Contains(y, S.Reals))
    elif isinstance(atom, (sp.GreaterThan, sp.LessThan)):
        x, y = atom.lhs, atom.rhs
        if y == 0:
            out.add(sp.Contains(x, S.Reals))
    elif isinstance(atom, sp.Contains):
        x, base = atom.args
        if base == S.Integers:
            out.update({sp.Contains(x, S.Rationals), sp.Contains(x, S.Reals)})
        elif base == S.Rationals:
            out.add(sp.Contains(x, S.Reals))
    return out


def close_assumptions(assumptions=True, *, max_rounds=3):
    original = normalize_assumptions(assumptions)
    facts = set(sp.And.make_args(original))
    for _ in range(max_rounds):
        expanded = set(facts)
        for fact in facts:
            expanded.update(_consequences(fact))
        if expanded == facts:
            break
        facts = expanded
    ordered = tuple(sorted(facts, key=sp.default_sort_key))
    return AssumptionClosure(
        original, frozenset(facts), sp.And(*ordered) if ordered else S.true
    )
