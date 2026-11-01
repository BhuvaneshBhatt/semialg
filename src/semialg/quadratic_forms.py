"""Exact constant-coefficient quadratic forms shared by structural rules."""

from dataclasses import dataclass

import sympy as sp

from .matrix_analysis import constant_symmetric_inertia


@dataclass(frozen=True)
class QuadraticForm:
    """The identity ``p(v) = v.T*A*v + b.T*v + c`` over real algebraics."""

    matrix: sp.ImmutableMatrix
    linear: sp.ImmutableMatrix
    constant: sp.Expr

    def stationary_point(self):
        """Return one stationary point, or None if the affine system is inconsistent."""
        try:
            solution, parameters = self.matrix.gauss_jordan_solve(-self.linear / 2, freevar=False)
        except ValueError:
            return None
        return solution.xreplace(dict.fromkeys(parameters, sp.S.Zero))


def quadratic_form(expression, variables):
    """Extract an exact real algebraic quadratic without differentiating expressions."""
    variables = tuple(variables)
    if not variables:
        return None
    try:
        poly = sp.Poly(expression, *variables)
    except sp.PolynomialError:
        return None
    if poly.total_degree() > 2 or any(
        coefficient.free_symbols
        or coefficient.is_real is not True
        or coefficient.is_algebraic is not True
        for coefficient in poly.coeffs()
    ):
        return None
    count = len(variables)
    matrix = sp.zeros(count)
    linear = sp.zeros(count, 1)
    constant = sp.S.Zero
    for powers, coefficient in poly.terms():
        support = [index for index, power in enumerate(powers) if power]
        if not support:
            constant = coefficient
        elif len(support) == 2:
            first, second = support
            matrix[first, second] = matrix[second, first] = coefficient / 2
        elif powers[support[0]] == 2:
            matrix[support[0], support[0]] = coefficient
        else:
            linear[support[0]] = coefficient
    return QuadraticForm(sp.ImmutableMatrix(matrix), sp.ImmutableMatrix(linear), constant)


def positive_quadratic_form(expression, variables):
    """Return a positive-definite quadratic form and its exact center and level."""
    form = quadratic_form(expression, variables)
    if form is None or constant_symmetric_inertia(form.matrix) != (len(variables), 0, 0):
        return None
    center = form.matrix.inv() * (-form.linear / 2)
    level = sp.cancel((center.T * form.matrix * center)[0] - form.constant)
    return form, center, level


def coercively_bounded(expression, variables):
    """Certify boundedness from coercive constraints and finite Boolean structure.

    A sublevel set of a coercive polynomial is bounded, including an empty
    sublevel. Subsets and finite unions retain boundedness. This routine only
    returns a positive certificate; failure to recognize a coercive constraint
    is not a proof that a region is unbounded.
    """
    if expression is sp.false or not variables:
        return True
    if isinstance(expression, sp.Or):
        return True if all(coercively_bounded(arg, variables) for arg in expression.args) else None
    if isinstance(expression, sp.And):
        return True if any(coercively_bounded(arg, variables) for arg in expression.args) else None
    if not isinstance(expression, sp.core.relational.Relational) or isinstance(
        expression, sp.Unequality
    ):
        return None
    residual = sp.expand(expression.lhs - expression.rhs)
    if isinstance(expression, (sp.GreaterThan, sp.StrictGreaterThan)):
        residual = -residual
    signs = (1, -1) if isinstance(expression, sp.Equality) else (1,)
    for sign in signs:
        try:
            poly = sp.Poly(sign * residual, *variables)
        except sp.PolynomialError:
            continue
        if any(c.free_symbols or c.is_real is not True for c in poly.coeffs()):
            continue
        degree = poly.total_degree()
        if degree <= 0 or degree % 2:
            continue
        if degree == 2:
            form = quadratic_form(poly.as_expr(), variables)
            if form is not None and constant_symmetric_inertia(form.matrix) == (
                len(variables),
                0,
                0,
            ):
                return True
            continue
        leading = sp.Add(
            *(
                coefficient * sp.prod(v**power for v, power in zip(variables, powers, strict=True))
                for powers, coefficient in poly.terms()
                if sum(powers) == degree
            )
        )
        norm = sp.Add(*(v * v for v in variables))
        ratio = sp.cancel(leading / (norm ** (degree // 2)))
        if not ratio.free_symbols and ratio.is_positive is True:
            return True
    return None
