from __future__ import annotations

import sympy as sp
from sympy import S

from ..exactness import certify_equal
from ..exceptions import POLYNOMIAL_ERRORS, SYMBOLIC_QUERY_ERRORS


def is_semialgebraic_expression(expr: sp.Expr, var: sp.Symbol) -> bool:
    """Conservative structural test for the subset handled by semialg.

    Algebraic powers with rational exponent, Abs/Min/Max/Piecewise and ordinary
    polynomial/rational arithmetic are accepted. Unknown Function heads are not.
    """
    expr = sp.sympify(expr)
    if not expr.has(var):
        return True
    if expr.is_Atom:
        return expr == var or expr.is_number or expr.is_Symbol
    if expr.func in (sp.Add, sp.Mul):
        return all(is_semialgebraic_expression(a, var) for a in expr.args)
    if expr.func is sp.Pow:
        b, e = expr.args
        return bool(e.is_Rational) and is_semialgebraic_expression(b, var)
    if expr.func in (sp.Abs, sp.Min, sp.Max):
        return all(is_semialgebraic_expression(a, var) for a in expr.args)
    if isinstance(expr, sp.Piecewise):
        return all(
            is_semialgebraic_expression(v, var) and is_semialgebraic_condition(c, var)
            for v, c in expr.args
        )
    return False


def is_semialgebraic_condition(cond, var):
    cond = sp.sympify(cond)
    if cond in (S.true, S.false):
        return True
    if isinstance(cond, (sp.And, sp.Or)):
        return all(is_semialgebraic_condition(a, var) for a in cond.args)
    if isinstance(cond, sp.Not):
        return is_semialgebraic_condition(cond.args[0], var)
    if isinstance(cond, sp.core.relational.Relational):
        return is_semialgebraic_expression(
            cond.lhs, var
        ) and is_semialgebraic_expression(cond.rhs, var)
    return False


def combine_domain(user_domain, natural_domain, var=None):
    """Intersect user and natural domains without losing the analysis variable.

    Set-valued domains require an explicit variable; using an unrelated dummy
    makes the restriction independent of the queried function variable.
    """
    if user_domain is None or user_domain is True or user_domain is S.true:
        return sp.sympify(natural_domain)
    if isinstance(user_domain, sp.Set):
        if var is None:
            raise ValueError("var is required when combining a Set-valued domain")
        return sp.And(sp.Contains(sp.sympify(var), user_domain), natural_domain)
    return sp.simplify(sp.And(sp.sympify(user_domain), sp.sympify(natural_domain)))


def condition_to_set(condition, var):
    """Convert ordinary univariate Boolean conditions to a SymPy Set.

    SymPy's ``Relational.as_set()`` can refuse parameterized endpoints such as
    ``x >= A`` even though the set is exactly ``[A, oo)`` with respect to x.
    Handle those one-variable affine endpoint relations structurally first.
    """
    condition = sp.sympify(condition)
    var = sp.sympify(var)
    if condition is S.true:
        return S.Reals
    if condition is S.false:
        return S.EmptySet
    if (
        isinstance(condition, sp.Contains)
        and condition.args[0] == var
        and isinstance(condition.args[1], sp.Set)
    ):
        return condition.args[1]
    if isinstance(condition, sp.And):
        pieces = [condition_to_set(a, var) for a in condition.args]
        if all(isinstance(q, sp.Set) for q in pieces):
            return sp.Intersection(*pieces)
    if isinstance(condition, sp.Or):
        pieces = [condition_to_set(a, var) for a in condition.args]
        if all(isinstance(q, sp.Set) for q in pieces):
            return sp.Union(*pieces)
    if isinstance(condition, sp.core.relational.Relational):
        lhs, rhs = condition.lhs, condition.rhs
        if lhs == var and not rhs.has(var):
            if isinstance(condition, sp.GreaterThan):
                return sp.Interval(rhs, sp.oo)
            if isinstance(condition, sp.StrictGreaterThan):
                return sp.Interval.open(rhs, sp.oo)
            if isinstance(condition, sp.LessThan):
                return sp.Interval(-sp.oo, rhs)
            if isinstance(condition, sp.StrictLessThan):
                return sp.Interval.open(-sp.oo, rhs)
            if isinstance(condition, sp.Equality):
                return sp.FiniteSet(rhs)
        if rhs == var and not lhs.has(var):
            if isinstance(condition, sp.GreaterThan):
                return sp.Interval(-sp.oo, lhs)
            if isinstance(condition, sp.StrictGreaterThan):
                return sp.Interval.open(-sp.oo, lhs)
            if isinstance(condition, sp.LessThan):
                return sp.Interval(lhs, sp.oo)
            if isinstance(condition, sp.StrictLessThan):
                return sp.Interval.open(lhs, sp.oo)
            if isinstance(condition, sp.Equality):
                return sp.FiniteSet(lhs)
    try:
        if isinstance(
            condition, (sp.logic.boolalg.Boolean, sp.core.relational.Relational)
        ):
            return condition.as_set()
    except (AttributeError, *SYMBOLIC_QUERY_ERRORS):
        # Semantic catalogue predicates are BooleanFunction subclasses without
        # SymPy's optional _eval_as_set implementation.
        pass
    return None


def set_to_condition(set_, y):
    """Convert a real Set to a Boolean membership condition in ``y``."""
    if set_ == S.Reals:
        return S.true
    if set_ == S.EmptySet:
        return S.false
    if isinstance(set_, sp.Interval):
        pieces = []
        if set_.start is not S.NegativeInfinity:
            pieces.append(y > set_.start if set_.left_open else y >= set_.start)
        if set_.end is not S.Infinity:
            pieces.append(y < set_.end if set_.right_open else y <= set_.end)
        return sp.And(*pieces) if pieces else S.true
    if isinstance(set_, sp.FiniteSet):
        return sp.Or(*(sp.Eq(y, v) for v in set_)) if set_ else S.false
    if isinstance(set_, sp.Union):
        return sp.Or(*(set_to_condition(a, y) for a in set_.args))
    if isinstance(set_, sp.Complement):
        a, b = set_.args
        return sp.And(set_to_condition(a, y), sp.Not(set_to_condition(b, y)))
    return sp.Contains(y, set_)


def interval_components(condition, var):
    s = condition_to_set(condition, var)
    if s is None:
        return None
    if isinstance(s, sp.Interval):
        return (s,)
    if isinstance(s, sp.Union) and all(isinstance(a, sp.Interval) for a in s.args):
        return tuple(s.args)
    if s == S.Reals:
        return (sp.Interval(-sp.oo, sp.oo),)
    return None


def is_full_natural_domain(requested_domain, natural_domain):
    if (
        requested_domain is True
        or requested_domain is S.true
        or requested_domain is None
    ):
        return True
    try:
        # Custom catalogue Boolean predicates must be lowered before generic
        # SymPy Boolean simplification; otherwise unknown predicate heads can
        # be simplified unsafely or merely syntactically.
        from ..stratification import semantic_condition

        req = semantic_condition(requested_domain)
        nat = semantic_condition(natural_domain)
        eq = sp.simplify(sp.simplify_logic(sp.Equivalent(req, nat)))
        return eq is S.true
    except SYMBOLIC_QUERY_ERRORS:
        return False


def exact_order_real_points(points):
    """Order exact real points only when every comparison is certified.

    Returns ``None``
    when a pair cannot be compared exactly.  This helper is shared by sign,
    monotonic-partition, extrema, and inverse reasoning.
    """
    out = []
    for q in points:
        q = sp.sympify(q)
        placed = False
        for i, r in enumerate(out):
            d = sp.simplify(q - r)
            if d.is_negative is True:
                out.insert(i, q)
                placed = True
                break
            if d.is_positive is True:
                continue
            if certify_equal(q, r) is True:
                placed = True
                break
            return None
        if not placed:
            out.append(q)
    return tuple(out)


def combine_domain_condition(natural_domain, user_domain, var):
    """Combine a natural-domain result with a requested real domain safely.

    A ``PropertyResult`` means the natural domain was not established as one
    Boolean condition.  In that case preserve the explicit user restriction,
    or use all reals when no restriction was supplied; never pass the result
    object itself into SymPy Boolean simplification.
    """
    var = sp.sympify(var)
    nat = None if hasattr(natural_domain, "primary") else sp.sympify(natural_domain)
    if user_domain is True or user_domain is None or user_domain is S.true:
        return S.true if nat is None else nat
    if isinstance(user_domain, sp.Set):
        usr = sp.Contains(var, user_domain)
    else:
        usr = sp.sympify(user_domain)
    if nat is None:
        return usr
    try:
        return sp.simplify(sp.And(nat, usr))
    except SYMBOLIC_QUERY_ERRORS:
        return sp.And(nat, usr)


def normalize_variables(var):
    """Return a nonempty tuple of distinct symbolic coordinates."""
    if isinstance(var, (tuple, list, sp.Tuple)):
        variables = tuple(sp.sympify(v) for v in var)
    else:
        variables = (sp.sympify(var),)
    if not variables:
        raise ValueError("at least one variable is required")
    if not all(isinstance(v, sp.Symbol) for v in variables):
        raise TypeError("variables must be Symbols")
    if len(set(variables)) != len(variables):
        raise ValueError("variables cannot contain duplicates")
    return variables


def normalize_expressions(expressions):
    """Return transformation components as a nonempty expression tuple."""
    if isinstance(expressions, (tuple, list, sp.Tuple, sp.MatrixBase)):
        values = tuple(sp.sympify(e) for e in expressions)
    else:
        values = (sp.sympify(expressions),)
    if not values:
        raise ValueError("at least one expression is required")
    return values


def normalize_domain_condition(domain, variables):
    """Return a Boolean domain condition for one or more coordinates."""
    variables = normalize_variables(variables)
    if domain in (None, True, S.true, S.Reals):
        return S.true
    if domain in (False, S.false, S.EmptySet):
        return S.false
    if isinstance(domain, sp.Set):
        if len(variables) == 1:
            return sp.Contains(variables[0], domain)
        return sp.Contains(sp.Tuple(*variables), domain)
    return sp.sympify(domain)


def normalize_domain_set(domain, variable):
    """Return a univariate real domain as a SymPy set when representable."""
    variable = normalize_variables(variable)[0]
    if domain in (None, True, S.true, S.Reals):
        return S.Reals
    if domain in (False, S.false, S.EmptySet):
        return S.EmptySet
    if isinstance(domain, sp.Set):
        return domain
    return condition_to_set(sp.sympify(domain), variable)


def affine_real_parameter_guard(expr, var):
    """Return ``(slope, offset, guard)`` for an affine real-valued expression.

    ``guard`` is parameter-only and is sufficient for ``expr`` to be real for
    every real ``var``.  ``None`` means that ``expr`` is not affine in ``var``.
    """
    expr, var = sp.sympify(expr), sp.sympify(var)
    try:
        p = sp.Poly(sp.expand(expr), var)
    except POLYNOMIAL_ERRORS:
        return None
    if p.degree() > 1:
        return None
    slope = sp.simplify(p.coeff_monomial(var))
    offset = sp.simplify(p.nth(0))
    real_slope = (
        S.true
        if slope.is_real is True
        else S.false
        if slope.is_real is False
        else sp.Q.real(slope)
    )
    real_offset = (
        S.true
        if offset.is_real is True
        else S.false
        if offset.is_real is False
        else sp.Q.real(offset)
    )
    return slope, offset, sp.And(real_slope, real_offset)
