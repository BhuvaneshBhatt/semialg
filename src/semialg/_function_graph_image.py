"""Exact graph/image QE backends for function-range computation."""

from __future__ import annotations

from collections.abc import Mapping

import sympy as sp
from sympy.logic.boolalg import BooleanFalse, BooleanTrue

from ._optimization_backends import (
    optimization_is_feasible as _optimization_is_feasible,
)
from ._optimization_backends import qe_by_complete_cad
from ._optimization_range import _EXPECTED_ERRORS
from ._range_special_cases import (
    _has_semialgebraic_special,
    _relation_for_function_graph,
)
from ._zero_testing import certified_nonzero, certified_zero
from .formula import parse_formula
from .function_graph import (
    UnsupportedFunctionGraph,
    semialgebraic_formula_graph,
    semialgebraic_function_graph,
)


def _free_image_coordinates(formula, bound_variables, value_symbol):
    # Coefficients outside the bound variables remain free image coordinates.
    parameters = formula.free_symbols - set(bound_variables) - {value_symbol}
    return (*sorted(parameters, key=sp.default_sort_key), value_symbol)


def _try_exact_graph_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Project a supported exact function graph onto its output coordinate."""

    try:
        graph = semialgebraic_function_graph(expr, value_symbol)
        condition_graph = semialgebraic_formula_graph(condition)
        image_formula = sp.And(condition_graph.formula, graph.formula)
        elimination_variables = tuple(
            dict.fromkeys(
                (
                    *variables,
                    *condition_graph.auxiliary_variables,
                    *graph.auxiliary_variables,
                )
            )
        )
        parsed = parse_formula(image_formula)
        free_variables = _free_image_coordinates(image_formula, elimination_variables, value_symbol)
        result = qe_by_complete_cad(
            (*free_variables, *elimination_variables),
            tuple(("exists", variable) for variable in elimination_variables),
            parsed,
            free_variables=free_variables,
            return_result=True,
        )
    except (UnsupportedFunctionGraph, *_EXPECTED_ERRORS):
        return None
    return sp.simplify(result.formula)


def _try_semialgebraic_graph_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Use CAD/QE on the shared graph for non-rational semialgebraic heads."""

    if not (_has_semialgebraic_special(expr) or _has_semialgebraic_special(condition)):
        return None
    return _try_exact_graph_image(expr, condition, variables, value_symbol)


def _substitute_formula(expr: sp.Expr, substitution: Mapping[sp.Symbol, sp.Expr]) -> sp.Expr:
    if expr is sp.true or isinstance(expr, BooleanTrue):
        return sp.true
    if expr is sp.false or isinstance(expr, BooleanFalse):
        return sp.false
    return sp.simplify(expr.subs(substitution))


def _try_affine_univariate_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variable: sp.Symbol,
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast exact image for an affine one-variable map.

    This covers important disconnected ranges without needing a full CAD run,
    e.g. the image of ``x`` over ``x <= -1 or x >= 1``.
    """

    try:
        poly = sp.Poly(sp.expand(expr), variable)
    except _EXPECTED_ERRORS:
        return None
    if poly.degree() != 1:
        return None
    coefficient = poly.coeff_monomial(variable)
    constant = poly.eval(0)
    if coefficient.is_real is False or constant.is_real is False:
        return None
    if certified_zero(coefficient) is True:
        sample = {variable: sp.Integer(0)}
        return (
            sp.Eq(value_symbol, sp.simplify(constant))
            if _optimization_is_feasible(condition, sample)
            else sp.false
        )
    # A coefficient that can vanish has a separate constant-map stratum.
    # Leave that family to parameter stratification or graph elimination.
    if not certified_nonzero(coefficient, assumptions=condition):
        return None
    inverse = sp.simplify((value_symbol - constant) / coefficient)
    return _substitute_formula(condition, {variable: inverse})


def _try_solved_graph_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Try exact image computation by solving the graph equation.

    Support fixed-coefficient univariate rational maps with rational inverse
    branches. Each branch retains its own poles before any substituted domain
    condition is simplified. Parameter-dependent degree loss and non-rational
    inverse branches require other image algorithms.
    """

    if len(variables) != 1:
        return None
    variable = variables[0]
    if not expr.is_rational_function(variable):
        return None
    affine = _try_affine_univariate_image(expr, condition, variable, value_symbol)
    if affine is not None:
        return sp.simplify(affine)
    if expr.free_symbols - {variable}:
        return None
    numerator, denominator = sp.fraction(sp.together(expr))
    equation = sp.expand(value_symbol * denominator - numerator)
    try:
        solutions = sp.solve(equation, variable)
    except _EXPECTED_ERRORS:
        return None
    if not solutions:
        return None
    branches: list[sp.Expr] = []
    for solution in solutions:
        if variable in solution.free_symbols or solution.has(sp.I):
            return None
        if any(
            isinstance(pow_expr, sp.Pow) and pow_expr.exp.is_integer is False
            for pow_expr in solution.atoms(sp.Pow)
        ):
            return None
        # For x = 1/y, Ne(x, 0) simplifies to True after substitution.
        # The inverse pole y = 0 must survive that simplification separately.
        _, inverse_denominator = sp.fraction(sp.together(solution))
        inverse_domain = sp.Ne(inverse_denominator, 0)
        branch_condition = _substitute_formula(condition, {variable: solution})
        branch_domain = _substitute_formula(sp.Ne(denominator, 0), {variable: solution})
        branches.append(sp.And(inverse_domain, branch_condition, branch_domain))
    return sp.simplify(sp.Or(*branches))


def _try_complete_cad_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Run the CAD/QE image formulation for polynomial/rational expressions."""

    relation, domain_constraint = _relation_for_function_graph(expr, value_symbol)
    image_formula = sp.And(condition, domain_constraint, relation)
    try:
        parsed = parse_formula(image_formula)
        free_variables = _free_image_coordinates(image_formula, variables, value_symbol)
        result = qe_by_complete_cad(
            (*free_variables, *variables),
            tuple(("exists", variable) for variable in variables),
            parsed,
            free_variables=free_variables,
            return_result=True,
        )
    except _EXPECTED_ERRORS:
        return None
    return sp.simplify(result.formula)


__all__ = []
