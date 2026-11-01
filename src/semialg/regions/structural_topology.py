"""Exact polynomial and Cartesian topology identities before CAD incidence."""

import sympy as sp

from ..exact_arithmetic import exact_sign, exact_truth
from ..matrix_analysis import constant_symmetric_inertia
from ..quadratic_forms import quadratic_form
from ..relations import split_relation


def _constant_topology(truth, operation):
    return sp.false if operation == "boundary" else sp.true if truth else sp.false


def polynomial_topology(atom, variables, operation):
    """Apply closed/open polynomial identities and complete quadratic sign geometry.

    For a quadratic, an indefinite Hessian or an inconsistent stationarity
    system makes every zero a sign-changing point. A semidefinite quadratic
    with a stationary point is classified by its unique critical value. This
    distinguishes an empty strict sublevel from a thin weak sublevel and a
    full weak superlevel, including singular quadratic forms.
    """
    if not isinstance(atom, sp.core.relational.Relational) or not variables:
        return None
    expression, operator = split_relation(atom)
    try:
        poly = sp.Poly(expression, *variables)
    except sp.PolynomialError:
        return None
    if any(c.free_symbols or c.is_real is not True for c in poly.coeffs()):
        return None
    if poly.total_degree() <= 0:
        return _constant_topology(exact_truth(atom), operation)
    zero = sp.Eq(atom.lhs, atom.rhs)
    if operator == "==":
        return sp.false if operation == "interior" else zero
    if operator == "!=":
        return sp.true if operation == "closure" else zero if operation == "boundary" else atom
    if operation == "closure" and operator in {"<=", ">="}:
        return atom
    if operation == "interior" and operator in {"<", ">"}:
        return atom
    form = quadratic_form(expression, variables)
    if form is None:
        return None
    inertia = constant_symmetric_inertia(form.matrix)
    if inertia is None:
        return None
    positive, negative, _zero = inertia
    center = None if positive and negative else form.stationary_point()
    if center is not None:
        critical = sp.cancel(form.constant + (form.linear.T * center)[0] / 2)
        # Reflect a nonpositive Hessian so the remaining classification has a
        # lower bound. No eigenvector or algebraic sample reconstruction is used.
        reflected_operator = operator
        if negative:
            critical = -critical
            reflected_operator = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}[operator]
        sign = exact_sign(critical)
        if sign > 0:
            return _constant_topology(reflected_operator in {">", ">="}, operation)
        if sign == 0:
            if reflected_operator == ">=":
                return _constant_topology(True, operation)
            if reflected_operator == "<":
                return sp.false
            if reflected_operator == "<=":
                return sp.false if operation == "interior" else zero
            return sp.true if operation == "closure" else zero if operation == "boundary" else atom
    if operation == "boundary":
        return zero
    if operation == "closure":
        constructor = sp.Le if operator == "<" else sp.Ge
        return constructor(atom.lhs, atom.rhs)
    constructor = sp.Lt if operator == "<=" else sp.Gt
    return constructor(atom.lhs, atom.rhs)


def cartesian_topology(expression, variables, operation):
    """Compute topology of a Cartesian product of exact degree-at-most-two 1-D sets.

    Each conjunction group involves a single coordinate, so exact univariate
    inequality reduction retains holes, equalities and strict endpoints. The
    closure and interior of a finite product are products of the corresponding
    factors; boundary is closure minus interior, also for thin or empty factors.
    """
    # The existing one-dimensional CAD path is cheaper than symbolic
    # inequality reduction for repeated interval-closure queries.
    if len(variables) < 2 or not isinstance(expression, sp.And):
        return None
    groups = {v: [] for v in variables}
    for atom in expression.args:
        if not isinstance(atom, sp.core.relational.Relational):
            return None
        support = atom.free_symbols
        if not support:
            if not exact_truth(atom):
                return sp.false
            continue
        if len(support) != 1 or not support.issubset(groups):
            return None
        variable = next(iter(support))
        polynomial = sp.Poly(atom.lhs - atom.rhs, variable)
        if polynomial.degree() > 2 or any(
            c.is_real is not True or c.is_algebraic is not True for c in polynomial.coeffs()
        ):
            return None
        groups[variable].append(atom)
    closures = []
    interiors = []
    boundaries = []
    for variable, atoms in groups.items():
        region = sp.reduce_inequalities(atoms, variable).as_set() if atoms else sp.S.Reals
        closures.append(region.closure.as_relational(variable))
        interiors.append(region.interior.as_relational(variable))
        boundaries.append(region.boundary.as_relational(variable))
    closure = sp.And(*closures)
    interior = sp.And(*interiors)
    if operation == "closure":
        return closure
    if operation == "interior":
        return interior
    return sp.Or(
        *(
            sp.And(boundary, *(item for j, item in enumerate(closures) if i != j))
            for i, boundary in enumerate(boundaries)
        )
    )


def structural_topology(expression, variables, operation):
    """Return an exact structural result, or None when CAD must decide."""
    if expression in (sp.true, sp.false):
        return _constant_topology(expression is sp.true, operation)
    try:
        result = polynomial_topology(expression, variables, operation)
        if result is not None:
            return result
        return cartesian_topology(expression, variables, operation)
    except (ValueError, TypeError, NotImplementedError, sp.PolynomialError):
        return None
