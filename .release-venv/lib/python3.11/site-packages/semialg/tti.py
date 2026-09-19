from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .ec.scoring import choose_eq_cons
from .formula import (
    Formula,
    Or,
    equational_constraints,
    formula_polynomials,
)


@dataclass(frozen=True)
class FormulaFamily:
    formula: Formula
    polynomials: tuple[sp.Expr, ...]
    equational_constraints: tuple[sp.Expr, ...]
    designated_ec: sp.Expr | None = None


def _flatten_or(formula: Formula) -> tuple[Formula, ...]:
    if isinstance(formula, Or):
        out: list[Formula] = []
        for arg in formula.args:
            out.extend(_flatten_or(arg))
        return tuple(out)
    return (formula,)


def extract_formula_families(formula: Formula) -> tuple[FormulaFamily, ...]:
    """
    Split a formula into TTICAD-style families.

    The strongest and most useful split for this implementation is across top-level OR
    branches, because truth-table invariance is typically exploited across lists of
    formulae or disjunctive branches. If there is no top-level OR, the entire matrix
    is treated as one family.
    """
    branches = _flatten_or(formula)
    families: list[FormulaFamily] = []
    for branch in branches:
        polys = tuple(dict.fromkeys(sp.expand(p) for p in formula_polynomials(branch)))
        ecs = tuple(dict.fromkeys(sp.expand(e) for e in equational_constraints(branch)))
        designated = choose_eq_cons(ecs, policy="lowest_degree")
        families.append(FormulaFamily(branch, polys, ecs, designated))
    return tuple(families)


def tti_eq_cons_by_level(
    formula: Formula, vars_: Sequence[sp.Symbol]
) -> dict[int, tuple[sp.Expr, ...]]:
    """
    A TTICAD-oriented EC map that keeps ECs branch-local, but merges them by level.
    """
    merged: dict[int, list[sp.Expr]] = {}
    for family in extract_formula_families(formula):
        for level in range(1, len(vars_) + 1):
            allowed = set(vars_[:level])
            current = [
                expr
                for expr in family.equational_constraints
                if expr.free_symbols and expr.free_symbols <= allowed
            ]
            if current:
                merged.setdefault(level, [])
                for expr in current:
                    if expr not in merged[level]:
                        merged[level].append(expr)
    return {level: tuple(exprs) for level, exprs in merged.items()}
