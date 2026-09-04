"""Fast exact special-case range and image reductions."""

from __future__ import annotations

from functools import cmp_to_key

import sympy as sp
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import BooleanFalse, BooleanTrue

from ._linear_relations import certified_sign, safe_linear_solution
from ._optimization_range import _EXPECTED_ERRORS, _finite_compare, _intervals_from_range_set
from .exact_arithmetic import compare_extended_reals
from .function_graph import has_semialgebraic_graph_special, semialgebraic_function_graph
from .relations import split_relation as _relation_parts


def _relation_for_function_graph(expr: sp.Expr, value_symbol: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    """Build the graph relation and denominator domain for a rational expression."""

    numerator, denominator = sp.fraction(sp.together(expr))
    numerator = sp.expand(numerator)
    denominator = sp.expand(denominator)
    relation = sp.Eq(sp.expand(value_symbol * denominator - numerator), 0)
    domain = sp.true if denominator == 1 else sp.Ne(denominator, 0)
    return relation, domain


def _is_supported_sqrt(expr: sp.Expr) -> bool:
    """Return whether ``expr`` is a direct principal square root."""

    return isinstance(expr, sp.Pow) and expr.exp == sp.Rational(1, 2)


def _is_semialgebraic_special(expr: sp.Expr) -> bool:
    """Return whether shared graph conversion needs a non-rational rule."""

    return has_semialgebraic_graph_special(expr)


def _has_semialgebraic_special(expr: sp.Expr) -> bool:
    return has_semialgebraic_graph_special(expr)


def _graph_formula_for_expression(
    expr: sp.Expr,
    target: sp.Symbol,
) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    """Return the shared exact real semialgebraic graph formula."""

    graph = semialgebraic_function_graph(expr, target)
    return graph.formula, graph.auxiliary_variables


def _range_formula_from_set(range_set: sp.Set, value_symbol: sp.Symbol) -> sp.Expr | None:
    """Convert a one-dimensional SymPy set into a range condition."""

    try:
        intervals = _intervals_from_range_set(range_set)
    except _EXPECTED_ERRORS:
        intervals = None
    if intervals is not None:
        pieces: list[sp.Expr] = []
        for left, right, left_closed, right_closed in intervals:
            if left == right:
                pieces.append(sp.Eq(value_symbol, left))
                continue
            lower = (
                sp.true
                if left == -sp.oo
                else (value_symbol >= left if left_closed else value_symbol > left)
            )
            upper = (
                sp.true
                if right == sp.oo
                else (value_symbol <= right if right_closed else value_symbol < right)
            )
            pieces.append(sp.And(lower, upper))
        return sp.Or(*pieces) if pieces else sp.false
    try:
        return range_set.as_relational(value_symbol)
    except _EXPECTED_ERRORS:
        return None


def _domain_set_from_condition(condition: sp.Expr, variable: sp.Symbol) -> sp.Set | None:
    """Return the one-variable domain represented by simple interval constraints."""

    if condition is sp.true or isinstance(condition, BooleanTrue):
        return sp.S.Reals
    if condition is sp.false or isinstance(condition, BooleanFalse):
        return sp.S.EmptySet
    atoms = condition.args if isinstance(condition, SymAnd) else (condition,)
    lower = -sp.oo
    upper = sp.oo
    left_open = False
    right_open = False
    for atom in atoms:
        if not getattr(atom, "is_Relational", False):
            return None
        expr, op = _relation_parts(atom)
        try:
            poly = sp.Poly(expr, variable)
        except _EXPECTED_ERRORS:
            return None
        if poly.degree() != 1:
            return None
        coeff = sp.simplify(poly.coeff_monomial(variable))
        const = sp.simplify(poly.eval(0))
        sign = certified_sign(coeff)
        if sign == 0:
            truth = bool(_relation_from_scalar(const, op))
            if not truth:
                return sp.S.EmptySet
            continue
        if op == "==":
            bound = safe_linear_solution(expr, variable)
            if bound is None:
                return None
        else:
            if sign not in (-1, 1):
                return None
            bound = safe_linear_solution(expr, variable)
            if bound is None:
                return None
            if sign < 0:
                op = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}.get(op, op)
        if op in {">", ">="}:
            if lower == -sp.oo or _finite_compare(lower, bound) < 0:
                lower = bound
                left_open = op == ">"
            elif sp.simplify(lower - bound) == 0:
                left_open = left_open or op == ">"
        elif op in {"<", "<="}:
            if upper == sp.oo or _finite_compare(bound, upper) < 0:
                upper = bound
                right_open = op == "<"
            elif sp.simplify(upper - bound) == 0:
                right_open = right_open or op == "<"
        elif op == "==":
            lower = upper = bound
            left_open = right_open = False
        else:
            return None
    if lower != -sp.oo and upper != sp.oo and _finite_compare(lower, upper) > 0:
        return sp.S.EmptySet
    return sp.Interval(lower, upper, left_open=left_open, right_open=right_open)


def _relation_from_scalar(value: sp.Expr, op: str) -> bool:
    if op == "<":
        return bool(value < 0)
    if op == "<=":
        return bool(value <= 0)
    if op == ">":
        return bool(value > 0)
    if op == ">=":
        return bool(value >= 0)
    if op == "==":
        return bool(sp.simplify(value) == 0)
    if op == "!=":
        return bool(sp.simplify(value) != 0)
    raise ValueError(op)


def _interval_from_bounds(lo: sp.Expr, hi: sp.Expr) -> sp.Interval:
    if _finite_compare(lo, hi) > 0:
        lo, hi = hi, lo
    return sp.Interval(sp.simplify(lo), sp.simplify(hi))


def _poly_range_on_interval(
    expr: sp.Expr, variable: sp.Symbol, domain: sp.Interval
) -> sp.Set | None:
    """Return the range of a univariate polynomial on a closed interval."""

    if not isinstance(domain, sp.Interval):
        return None
    try:
        poly = sp.Poly(sp.expand(expr), variable)
    except _EXPECTED_ERRORS:
        return None
    if domain.start in (-sp.oo, sp.oo) or domain.end in (-sp.oo, sp.oo):
        return None
    candidates = [domain.start, domain.end]
    derivative = sp.diff(poly.as_expr(), variable)
    try:
        for root in sp.solve(sp.Eq(derivative, 0), variable):
            if root.is_real is False:
                continue
            if bool(root >= domain.start) and bool(root <= domain.end):
                candidates.append(root)
    except _EXPECTED_ERRORS:
        pass
    values = [sp.simplify(expr.subs(variable, item)) for item in candidates]
    lo = min(values, key=cmp_to_key(_compare_bounds))
    hi = max(values, key=cmp_to_key(_compare_bounds))
    return _interval_from_bounds(lo, hi)


def _range_via_abs_affine(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for affine expressions in ``Abs(variable)``."""

    if len(variables) != 1:
        return None
    variable = variables[0]
    abs_expr = sp.Abs(variable)
    if not expr.has(abs_expr):
        return None
    domain = _domain_set_from_condition(condition, variable)
    if domain is None:
        return None
    u = sp.Symbol("_abs_value", real=True)
    transformed = expr.xreplace({abs_expr: u})
    try:
        poly = sp.Poly(sp.expand(transformed), u)
    except _EXPECTED_ERRORS:
        return None
    if poly.degree() != 1 or any(sym != u for sym in transformed.free_symbols):
        return None
    coeff = sp.simplify(poly.coeff_monomial(u))
    const = sp.simplify(poly.eval(0))
    if coeff == 0:
        return _range_formula_from_set(sp.FiniteSet(const), value_symbol)
    if domain is sp.S.Reals:
        abs_range = sp.Interval(0, sp.oo)
    elif isinstance(domain, sp.Interval):
        values = []
        if domain.start != -sp.oo:
            values.append(sp.Abs(domain.start))
        if domain.end != sp.oo:
            values.append(sp.Abs(domain.end))
        crosses_zero = domain.start <= 0 and domain.end >= 0
        lower = (
            sp.Integer(0) if bool(crosses_zero) else min(values, key=cmp_to_key(_compare_bounds))
        )
        upper = (
            sp.oo
            if domain.start == -sp.oo or domain.end == sp.oo
            else max(values, key=cmp_to_key(_compare_bounds))
        )
        abs_range = sp.Interval(lower, upper)
    else:
        return None
    left = sp.simplify(coeff * abs_range.start + const) if abs_range.start != -sp.oo else -sp.oo
    if abs_range.end != sp.oo:
        right = sp.simplify(coeff * abs_range.end + const)
    else:
        sign = certified_sign(coeff)
        if sign == 1:
            right = sp.oo
        elif sign == -1:
            right = -sp.oo
        else:
            return None
    return _range_formula_from_set(_interval_from_bounds(left, right), value_symbol)


def _range_via_sqrt_quadratic(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for square roots of nonnegative quadratic radicands."""

    if len(variables) != 1 or not _is_supported_sqrt(expr):
        return None
    variable = variables[0]
    domain = _domain_set_from_condition(condition, variable)
    if domain is None:
        return None
    try:
        base_poly = sp.Poly(sp.expand(expr.base), variable)
    except _EXPECTED_ERRORS:
        return None
    if base_poly.degree() > 2:
        return None
    coeff2 = sp.simplify(base_poly.coeff_monomial(variable**2))
    coeff1 = sp.simplify(base_poly.coeff_monomial(variable))
    coeff0 = sp.simplify(base_poly.coeff_monomial(1))
    if coeff2 == -1 and coeff1 == 0 and coeff0.is_real:
        natural = sp.Interval(-sp.sqrt(coeff0), sp.sqrt(coeff0))
    else:
        return None
    active_domain = natural if domain is sp.S.Reals else domain.intersect(natural)
    if not isinstance(active_domain, sp.Interval):
        return None
    base_range = _poly_range_on_interval(expr.base, variable, active_domain)
    if not isinstance(base_range, sp.Interval):
        return None
    upper = sp.sqrt(sp.simplify(base_range.end))
    return _range_formula_from_set(sp.Interval(0, upper), value_symbol)


def _range_via_sympy_calculus(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast image computation for common univariate expression heads."""

    if len(variables) != 1:
        return None
    variable = variables[0]
    domain = _domain_set_from_condition(condition, variable)
    if domain is None:
        return None
    try:
        from sympy.calculus.util import function_range as sympy_range

        range_set = sympy_range(expr, variable, domain)
    except _EXPECTED_ERRORS:
        return None
    return _range_formula_from_set(range_set, value_symbol)


def _compare_bounds(left: sp.Expr, right: sp.Expr) -> int:
    """Compare extended-real range bounds without numerical approximation."""

    return compare_extended_reals(left, right)


def _range_via_max_min(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for two-argument Max/Min with affine arguments."""

    if len(variables) != 1 or expr.func not in {sp.Max, sp.Min} or len(expr.args) != 2:
        return None
    variable = variables[0]
    domain = _domain_set_from_condition(condition, variable)
    if not isinstance(domain, sp.Interval):
        return None
    arg_ranges: list[sp.Set] = []
    for arg in expr.args:
        try:
            poly = sp.Poly(sp.expand(arg), variable)
        except _EXPECTED_ERRORS:
            return None
        if poly.degree() > 1:
            return None
        coeff = sp.simplify(poly.coeff_monomial(variable))
        const = sp.simplify(poly.eval(0))
        if coeff == 0:
            arg_ranges.append(sp.FiniteSet(const))
            continue
        left = -sp.oo if domain.start == -sp.oo else sp.simplify(coeff * domain.start + const)
        right = sp.oo if domain.end == sp.oo else sp.simplify(coeff * domain.end + const)
        lo, hi = (left, right) if _finite_compare(left, right) <= 0 else (right, left)
        arg_ranges.append(sp.Interval(lo, hi))
    if expr.func is sp.Max:
        lower = max((item.inf for item in arg_ranges), key=cmp_to_key(_compare_bounds))
        upper = max((item.sup for item in arg_ranges), key=cmp_to_key(_compare_bounds))
    else:
        lower = min((item.inf for item in arg_ranges), key=cmp_to_key(_compare_bounds))
        upper = min((item.sup for item in arg_ranges), key=cmp_to_key(_compare_bounds))
    return _range_formula_from_set(sp.Interval(lower, upper), value_symbol)


def _range_via_piecewise(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for univariate Piecewise expressions by branch images."""

    if len(variables) != 1 or not isinstance(expr, sp.Piecewise):
        return None
    variable = variables[0]
    pieces: list[sp.Set] = []
    previous: list[sp.Expr] = []
    for branch_expr, branch_cond in expr.args:
        effective = sp.And(condition, branch_cond, *(sp.Not(item) for item in previous))
        previous.append(branch_cond)
        domain = _domain_set_from_condition(effective, variable)
        if domain is None:
            try:
                reduced = sp.reduce_inequalities(
                    list(effective.args) if isinstance(effective, SymAnd) else [effective],
                    variable,
                )
                domain = reduced.as_set()
            except _EXPECTED_ERRORS:
                return None
        if domain is sp.S.EmptySet:
            continue
        if isinstance(domain, sp.Union):
            domains = domain.args
        else:
            domains = (domain,)
        for branch_domain in domains:
            image = _poly_range_on_interval(branch_expr, variable, branch_domain)
            if image is None:
                return None
            pieces.append(image)
    if not pieces:
        return sp.false
    return _range_formula_from_set(sp.Union(*pieces), value_symbol)


def _try_direct_special_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Use fast univariate range formulas before constructing a CAD image problem."""

    for helper in (
        _range_via_abs_affine,
        _range_via_sqrt_quadratic,
        _range_via_max_min,
        _range_via_piecewise,
        _range_via_sympy_calculus,
    ):
        result = helper(expr, condition, variables, value_symbol)
        if result is not None:
            return result
    return None


__all__ = []
