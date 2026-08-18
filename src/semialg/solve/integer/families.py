from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp
from sympy import Eq


@dataclass(frozen=True)
class IntegerFamilyTag:
    name: str
    metadata: dict


def detect_sum_fam(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> IntegerFamilyTag | None:
    variables = tuple(variables)
    if len(variables) != 2:
        return None
    atoms = list(expr.args) if isinstance(expr, sp.And) else [expr]
    eqs = [a for a in atoms if isinstance(a, Eq)]
    others = [a for a in atoms if not isinstance(a, Eq)]
    if len(eqs) != 1 or others:
        return None
    x, y = variables
    poly = sp.Poly(sp.expand(eqs[0].lhs - eqs[0].rhs), x, y)
    if (
        poly.coeff_monomial(x**2) == 1
        and poly.coeff_monomial(y**2) == 1
        and poly.coeff_monomial(x * y) == 0
    ):
        c = -poly.coeff_monomial(1)
        if c.is_integer:
            return IntegerFamilyTag("sum_of_two_squares", {"target_norm": int(c)})
    return None


def detect_binary_homog_fam(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> IntegerFamilyTag | None:
    variables = tuple(variables)
    if len(variables) != 2:
        return None
    atoms = list(expr.args) if isinstance(expr, sp.And) else [expr]
    eqs = [a for a in atoms if isinstance(a, Eq)]
    others = [a for a in atoms if not isinstance(a, Eq)]
    if len(eqs) != 1 or others:
        return None
    x, y = variables
    poly = sp.Poly(sp.expand(eqs[0].lhs - eqs[0].rhs), x, y)
    degs = {sum(mon) for mon, coeff in poly.terms() if coeff != 0}
    if len(degs) == 1:
        return IntegerFamilyTag("binary_homogeneous", {"degree": list(degs)[0], "polynomial": poly})
    return None


def detect_int_problem_fam(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> IntegerFamilyTag | None:
    for detector in (
        detect_sum_fam,
        detect_binary_homog_fam,
    ):
        tag = detector(expr, variables)
        if tag is not None:
            return tag
    return None


__all__ = [
    "IntegerFamilyTag",
    "detect_sum_fam",
    "detect_binary_homog_fam",
    "detect_int_problem_fam",
]
