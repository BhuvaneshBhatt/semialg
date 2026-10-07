"""Bounded parameter-cell refinement for conditional results."""

from __future__ import annotations

import sympy as sp
from sympy import S

from .assumption_closure import close_assumptions
from .exactness import certify_identity_equal
from .results import ConditionalExpression, PropertyResult
from .stratification import normalize_cases, semantic_condition


def minimize_strata(result):
    if not isinstance(result, PropertyResult) or not result.cases:
        return result
    cases = []
    for case in result.cases:
        cond = semantic_condition(case.condition)
        atoms = list(sp.And.make_args(cond))
        kept = []
        for i, atom in enumerate(atoms):
            others = sp.And(*(atoms[:i] + atoms[i + 1 :])) if len(atoms) > 1 else S.true
            if close_assumptions(others).entails(atom) is True:
                continue
            kept.append(atom)
        cases.append(
            ConditionalExpression(case.value, sp.And(*kept) if kept else S.true)
        )
    return normalize_cases(
        cases,
        source=result.source,
        primary=result.primary,
        all_properties=result.all_properties,
        evidence=result.evidence,
    )


def specialize_cells(result, substitutions):
    if not isinstance(result, PropertyResult):
        return result
    cases = []
    for c in result.cases:
        cond = sp.simplify(c.condition.subs(substitutions))
        if cond is not S.false:
            cases.append(ConditionalExpression(c.value, cond))
    return normalize_cases(
        cases, source=result.source, primary=result.primary, evidence=result.evidence
    )


def refine_parameter_cells(
    result, parameters=None, *, split_zero=True, split_sign=False, max_cells=32
):
    """Refine unresolved parameter cells by exact zero/sign partitions.

    This is deliberately bounded and purely symbolic.  It never samples parameter
    values and never labels a cell with a new conclusion; it only refines guards so
    existing conditional conclusions can specialize cleanly downstream.
    """
    if not isinstance(result, PropertyResult) or not result.cases:
        return result
    if parameters is None:
        parameters = sorted(
            set().union(*(c.condition.free_symbols for c in result.cases)),
            key=sp.default_sort_key,
        )
    guards = [S.true]
    for p in parameters:
        pieces = [sp.Eq(p, 0), sp.Ne(p, 0)] if split_zero else [S.true]
        if split_sign:
            pieces = [p < 0, sp.Eq(p, 0), p > 0]
        new = []
        for g in guards:
            for q in pieces:
                h = sp.simplify_logic(sp.And(g, q))
                if h is not S.false:
                    new.append(h)
                if len(new) >= max_cells:
                    break
            if len(new) >= max_cells:
                break
        guards = new or guards
        if len(guards) >= max_cells:
            break
    out = []
    for case in result.cases:
        for g in guards:
            h = sp.simplify_logic(sp.And(case.condition, g))
            if h is not S.false:
                out.append(ConditionalExpression(case.value, h))
    return normalize_cases(
        out,
        source=result.source,
        primary=result.primary,
        all_properties=result.all_properties,
        evidence=result.evidence,
    )


def exponential_affine_transition_cells(expr, variable, parameter):
    """Certified simple parameter strata for ``exp(x) + a``-type translations.

    This intentionally recognizes only a family whose transition is exact: exp(x)+a
    has one real root for a<0 and none for a>=0 (with a=0 approached only at -oo).
    """
    x = sp.sympify(variable)
    a = sp.sympify(parameter)
    f = sp.expand(sp.sympify(expr))
    if certify_identity_equal(f, sp.exp(x) + a) is not True:
        return None
    return PropertyResult(
        None,
        cases=(ConditionalExpression(1, a < 0), ConditionalExpression(0, a >= 0)),
        source="exponential-translation-stratification",
    )


def quadratic_discriminant_cells(expr, variable, parameter):
    """Exact parameter root-count strata for univariate quadratics in x."""
    x = sp.sympify(variable)
    p = sp.sympify(parameter)
    try:
        poly = sp.Poly(expr, x)
    except sp.PolynomialError:
        return None
    if poly.degree() != 2 or p not in poly.as_expr().free_symbols:
        return None
    A, B, C = poly.all_coeffs()
    D = sp.factor(B**2 - 4 * A * C)
    return PropertyResult(
        None,
        cases=(
            ConditionalExpression(2, D > 0),
            ConditionalExpression(1, sp.Eq(D, 0)),
            ConditionalExpression(0, D < 0),
        ),
        source="quadratic-discriminant-stratification",
        evidence=(D,),
    )
