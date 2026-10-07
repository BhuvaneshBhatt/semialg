"""Normalization of conditional inputs for public property queries."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from .assumption_normalization import normalize_assumptions
from .real_analysis.common import normalize_domain_condition, normalize_variables
from .results import ConditionalExpression


@dataclass(frozen=True)
class NormalizedConditionalInput:
    expr: object
    condition: sp.Basic = S.true


def normalize_conditional_input(expr) -> NormalizedConditionalInput:
    """Flatten nested :class:`ConditionalExpression` inputs without symbol capture."""
    conditions = []
    value = expr
    while isinstance(value, ConditionalExpression):
        conditions.append(sp.sympify(value.condition))
        value = value.value
    condition = sp.And(*conditions) if conditions else S.true
    return NormalizedConditionalInput(value, sp.simplify(condition))


def restrict_domain(domain, condition, var):
    """Intersect a query domain with a condition without changing dimension."""
    condition = sp.sympify(condition)
    if condition is S.true:
        return domain
    variables = normalize_variables(var)
    base = normalize_domain_condition(domain, variables)
    if base is S.true:
        return condition
    if base is S.false:
        return S.false
    return sp.And(base, condition)


def combine_assumptions(assumptions, condition):
    condition = sp.sympify(condition)
    if condition is S.true:
        return assumptions
    return normalize_assumptions(sp.And(sp.sympify(assumptions), condition))
