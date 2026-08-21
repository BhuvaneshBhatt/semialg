from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .formula_utils import conjuncts as _conjuncts


@dataclass(frozen=True)
class LinDivisReduction:
    solved_variable: sp.Symbol
    denominator: sp.Expr
    numerator: sp.Expr
    replacement: sp.Expr
    divisibility_condition: sp.Expr
    reduced_formula: sp.Expr


def detect_lin_reduction(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> LinDivisReduction | None:
    variables = tuple(variables)
    conjuncts = _conjuncts(expr)
    for atom in conjuncts:
        if not isinstance(atom, sp.Equality):
            continue
        diff = sp.expand(atom.lhs - atom.rhs)
        for var in reversed(variables):
            poly = sp.Poly(diff, var)
            if poly.degree() != 1:
                continue
            coeff = sp.expand(poly.coeff_monomial(var))
            const = sp.expand(poly.coeff_monomial(1))
            if coeff == 0 or coeff.has(var) or const.has(var):
                continue
            numerator = sp.expand(-const)
            denominator = sp.expand(coeff)
            replacement = sp.simplify(numerator / denominator)
            divisibility = sp.Eq(sp.Mod(numerator, denominator), 0)
            rest = [a for a in conjuncts if a is not atom]
            substituted = [sp.simplify(a.subs(var, replacement)) for a in rest]
            reduced_formula = sp.And(divisibility, *substituted) if substituted else divisibility
            return LinDivisReduction(
                solved_variable=var,
                denominator=denominator,
                numerator=numerator,
                replacement=replacement,
                divisibility_condition=divisibility,
                reduced_formula=sp.simplify(reduced_formula),
            )
    return None


def apply_lin_reduction(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, LinDivisReduction | None]:
    reduction = detect_lin_reduction(expr, variables)
    if reduction is None:
        return expr, None
    reconstructed = sp.And(
        sp.Eq(reduction.solved_variable, reduction.replacement),
        reduction.reduced_formula,
    )
    return sp.simplify(reconstructed), reduction


__all__ = [
    "LinDivisReduction",
    "detect_lin_reduction",
    "apply_lin_reduction",
]
