"""Bounded Boolean Fourier–Motzkin elimination before polynomial backends."""

from dataclasses import dataclass
from itertools import product

import sympy as sp
from sympy.core.relational import Relational

from ..presolve import fourier_motzkin_eliminate, presolve_semialgebraic


@dataclass(frozen=True)
class LinearEliminationResult:
    """Exact linear elimination with a bounded Boolean expansion."""

    formula: sp.Expr
    method: str = "linear-boolean"


def _conjunctions(expression, limit):
    """Expand NNF with a cap checked before taking each Cartesian product."""
    if expression is sp.true:
        return [()]
    if expression is sp.false:
        return []
    if isinstance(expression, sp.Or):
        result = []
        for child in expression.args:
            terms = _conjunctions(child, limit)
            if terms is None or len(result) + len(terms) > limit:
                return None
            result.extend(terms)
        return result
    if isinstance(expression, sp.And):
        result = [()]
        for child in expression.args:
            terms = _conjunctions(child, limit)
            if terms is None or len(result) * len(terms) > limit:
                return None
            result = [left + right for left, right in product(result, terms)]
        return result
    if isinstance(expression, sp.Equality):
        return [(expression.lhs <= expression.rhs, expression.lhs >= expression.rhs)]
    if isinstance(expression, sp.Unequality):
        if limit < 2:
            return None
        return [(expression.lhs < expression.rhs,), (expression.lhs > expression.rhs,)]
    if isinstance(expression, Relational):
        return [(expression,)]
    return None


def try_linear_elimination(matrix, quantifiers, variables, *, max_terms=64):
    """Eliminate an affine Boolean formula, including alternating prefixes.

    Universal blocks use negation duality; existential quantification
    distributes over the bounded disjunction. Coefficients must be exact real
    constants, so parameter-dependent leading coefficients retain their zero
    strata in the general solver. Declining leaves the original problem intact.
    """
    if not quantifiers:
        return None
    for atom in matrix.atoms(Relational):
        try:
            poly = sp.Poly(atom.lhs - atom.rhs, *variables, extension=True)
        except (sp.PolynomialError, ValueError, TypeError):
            return None
        if poly.total_degree() > 1:
            return None
        if not all(coefficient.is_real is True for coefficient in poly.coeffs()):
            return None
    current = matrix
    for kind, variable in reversed(tuple(quantifiers)):
        universal = kind.lower() == "forall"
        if kind.lower() not in {"exists", "forall"}:
            raise ValueError(f"unsupported quantifier: {kind!r}")
        expression = sp.to_nnf(sp.Not(current) if universal else current)
        terms = _conjunctions(expression, max_terms)
        if terms is None:
            return None
        images = []
        for term in terms:
            conjunction = sp.And(*term)
            if variable not in conjunction.free_symbols:
                images.append(conjunction)
                continue
            presolved = presolve_semialgebraic(conjunction, variables, eliminate=(variable,))
            if variable not in presolved.formula.free_symbols:
                images.append(presolved.formula)
                continue
            image = fourier_motzkin_eliminate(presolved.formula, (variable,))
            if image is None:
                return None
            images.append(image[0])
        current = sp.Or(*images)
        if universal:
            current = sp.to_nnf(sp.Not(current))
    return LinearEliminationResult(current)
