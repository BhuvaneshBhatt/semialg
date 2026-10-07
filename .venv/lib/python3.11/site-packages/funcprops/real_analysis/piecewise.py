"""Certified decomposition and merging support for real ``Piecewise`` expressions.

A Piecewise expression is an ordered decision list: later conditions only apply
where all earlier conditions are false.  This module turns that semantics into
explicit effective branch cells, records which branch owns every finite seam,
and certifies one-sided value/derivative compatibility without asking SymPy for
a limit of the whole Piecewise expression.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..exactness import certify_equal
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import Boundedness, Monotonicity, Sign
from .common import condition_to_set


@dataclass(frozen=True)
class PiecewiseBranch:
    """One connected effective branch cell.

    ``source_index`` identifies the original Piecewise arm.  ``condition`` is
    the normalized *effective* condition after excluding all earlier arms.
    ``domain`` is the corresponding real set/cell when it can be represented.
    """

    expression: sp.Expr
    condition: sp.Basic
    domain: sp.Set | None
    source_index: int = 0
    source_condition: sp.Basic = S.true


@dataclass(frozen=True)
class PiecewiseBoundary:
    """A finite seam between two effective branch cells."""

    point: sp.Expr
    left_branch: int | None
    right_branch: int | None
    owner_branch: int | None


@dataclass(frozen=True)
class PiecewiseBoundaryCertificate:
    """Certified local data at one Piecewise seam."""

    boundary: PiecewiseBoundary
    value: sp.Expr | None
    left_limit: sp.Expr | None
    right_limit: sp.Expr | None
    derivative_limits: tuple[tuple[sp.Expr | None, sp.Expr | None], ...] = ()
    continuous: bool | None = None
    differentiability_order: int | None = None
    complete: bool = False


@dataclass(frozen=True)
class PiecewiseDecomposition:
    """Effective ordered branch cells and their finite boundaries."""

    expression: sp.Piecewise
    variable: sp.Symbol
    requested_domain: sp.Basic | sp.Set
    branches: tuple[PiecewiseBranch, ...]
    boundaries: tuple[PiecewiseBoundary, ...]
    complete: bool


def _effective_condition(domain, condition, covered, var):
    domain_condition = (
        domain.as_relational(var) if isinstance(domain, sp.Set) else sp.sympify(domain)
    )
    try:
        return sp.simplify_logic(sp.And(domain_condition, condition, sp.Not(covered)))
    except SYMBOLIC_QUERY_ERRORS:
        return sp.And(domain_condition, condition, sp.Not(covered))


def _components(domain_set):
    if isinstance(domain_set, sp.Union):
        return tuple(domain_set.args)
    return (domain_set,)


def _condition_for_component(component, var):
    if isinstance(component, sp.Interval):
        return component.as_relational(var)
    if isinstance(component, sp.FiniteSet):
        vals = tuple(component)
        return sp.Or(*(sp.Eq(var, v) for v in vals)) if vals else S.false
    return sp.Contains(var, component)


def piecewise_branches(expr, var, domain=S.true):
    """Return nonempty connected effective branches of a ``Piecewise`` expression."""
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if not isinstance(expr, sp.Piecewise):
        return ()
    covered = S.false
    branches = []
    for index, (value, condition) in enumerate(expr.args):
        condition = sp.sympify(condition)
        effective = _effective_condition(domain, condition, covered, var)
        covered = sp.Or(covered, condition)
        if effective is S.false:
            continue
        domain_set = condition_to_set(effective, var)
        if domain_set is S.EmptySet:
            continue
        if domain_set is None:
            branches.append(
                PiecewiseBranch(sp.sympify(value), effective, None, index, condition)
            )
            continue
        for component in _components(domain_set):
            if component is S.EmptySet:
                continue
            component_condition = _condition_for_component(component, var)
            branches.append(
                PiecewiseBranch(
                    sp.sympify(value),
                    sp.simplify(sp.And(effective, component_condition)),
                    component,
                    index,
                    condition,
                )
            )
    return tuple(branches)


def _contains_point(branch, point, var):
    if branch.domain is not None:
        try:
            hit = branch.domain.contains(point)
            if hit in (S.true, True):
                return True
            if hit in (S.false, False):
                return False
        except SYMBOLIC_QUERY_ERRORS:
            pass
    try:
        verdict = sp.simplify(branch.condition.subs(var, point))
    except SYMBOLIC_QUERY_ERRORS:
        return None
    if verdict in (S.true, True):
        return True
    if verdict in (S.false, False):
        return False
    return None


def _interval_branches(branches):
    return tuple(b for b in branches if isinstance(b.domain, sp.Interval))


def _branch_index(branches, target):
    for i, branch in enumerate(branches):
        if branch is target:
            return i
    return None


def _extract_boundaries(branches, var):
    interval_branches = sorted(
        _interval_branches(branches), key=lambda b: sp.default_sort_key(b.domain.start)
    )
    points = set()
    for b in interval_branches:
        if b.domain.start not in (S.NegativeInfinity, -S.Infinity):
            points.add(b.domain.start)
        if b.domain.end not in (S.Infinity, S.NegativeInfinity):
            points.add(b.domain.end)
    out = []
    for point in sorted(points, key=sp.default_sort_key):
        left = None
        right = None
        for b in interval_branches:
            if certify_equal(b.domain.end, point) is True:
                left = _branch_index(branches, b)
            if certify_equal(b.domain.start, point) is True and right is None:
                right = _branch_index(branches, b)
        # A boundary at only the outside edge of the requested domain is not an
        # interface unless some isolated/other cell owns that point.
        candidates = [
            i for i, b in enumerate(branches) if _contains_point(b, point, var) is True
        ]
        owner = (
            min(candidates, key=lambda i: branches[i].source_index)
            if candidates
            else None
        )
        if left is not None and right is not None and left != right:
            out.append(PiecewiseBoundary(point, left, right, owner))
        elif owner is not None and (
            (left is not None and owner != left)
            or (right is not None and owner != right)
        ):
            out.append(PiecewiseBoundary(point, left, right, owner))
    return tuple(out)


def piecewise_decomposition(expr, var, domain=S.true):
    """Build the certified effective-branch decomposition for a real query."""
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if not isinstance(expr, sp.Piecewise):
        raise TypeError("piecewise_decomposition requires a SymPy Piecewise expression")
    branches = piecewise_branches(expr, var, domain)
    complete = bool(branches) and all(b.domain is not None for b in branches)
    boundaries = _extract_boundaries(branches, var) if complete else ()
    return PiecewiseDecomposition(expr, var, domain, branches, boundaries, complete)


def _side_limit(branch, var, point, direction, derivative_order=0):
    from ..proof_cache import cached_derivative, cached_endpoint_limit

    if branch is None:
        return None
    value = branch.expression
    if derivative_order:
        try:
            value = cached_derivative(value, var, order=derivative_order)
        except SYMBOLIC_QUERY_ERRORS:
            return None
    return cached_endpoint_limit(value, var, point, direction=direction, attained=False)


def certify_piecewise_boundaries(decomposition, *, max_order=1):
    """Certify seam values and one-sided derivatives through ``max_order``.

    Limits are always taken on the selected branch expression, never on the
    enclosing Piecewise object.  ``differentiability_order`` is the highest
    derivative order whose two one-sided limits agree at the seam (0 means
    continuity only).
    """
    if not decomposition.complete:
        return ()
    branches = decomposition.branches
    var = decomposition.variable
    out = []
    for boundary in decomposition.boundaries:
        left = (
            branches[boundary.left_branch] if boundary.left_branch is not None else None
        )
        right = (
            branches[boundary.right_branch]
            if boundary.right_branch is not None
            else None
        )
        owner = (
            branches[boundary.owner_branch]
            if boundary.owner_branch is not None
            else None
        )
        try:
            value = (
                owner.expression.subs(var, boundary.point)
                if owner is not None
                else decomposition.expression.subs(var, boundary.point)
            )
        except SYMBOLIC_QUERY_ERRORS:
            value = None
        left0 = _side_limit(left, var, boundary.point, "-")
        right0 = _side_limit(right, var, boundary.point, "+")
        continuous = None
        if value is not None and left0 is not None and right0 is not None:
            lv = certify_equal(left0, value)
            rv = certify_equal(right0, value)
            if lv is True and rv is True:
                continuous = True
            elif lv is False or rv is False:
                continuous = False
        derivative_limits = []
        order = 0 if continuous is True else None
        if continuous is True:
            for k in range(1, max_order + 1):
                lval = _side_limit(left, var, boundary.point, "-", k)
                rval = _side_limit(right, var, boundary.point, "+", k)
                derivative_limits.append((lval, rval))
                if lval is None or rval is None:
                    break
                delta = certify_equal(lval, rval)
                if delta is True:
                    order = k
                    continue
                if delta is False:
                    break
                order = None
                break
        complete = value is not None and left0 is not None and right0 is not None
        out.append(
            PiecewiseBoundaryCertificate(
                boundary,
                value,
                left0,
                right0,
                tuple(derivative_limits),
                continuous,
                order,
                complete,
            )
        )
    return tuple(out)


def combine_piecewise_sets(values):
    """Merge branch-local exact sets without losing openness/endpoints."""
    all_values = tuple(values)
    vals = tuple(v for v in all_values if isinstance(v, sp.Set))
    if len(vals) != len(all_values):
        return None
    return sp.Union(*vals) if vals else S.EmptySet


def combine_piecewise_conditions(values, variable=None):
    """Merge branch-local image conditions, using exact set union when possible."""
    vals = tuple(sp.sympify(v) for v in values)
    if not vals:
        return S.false
    if variable is None:
        symbols = set().union(*(v.free_symbols for v in vals))
        variable = next(iter(symbols)) if len(symbols) == 1 else None
    if variable is not None:
        sets = tuple(condition_to_set(v, variable) for v in vals)
        if all(isinstance(q, sp.Set) for q in sets):
            union = sp.Union(*sets)
            if union == S.Reals:
                return S.true
            if union == S.EmptySet:
                return S.false
            # Preserve branch-local witnesses for ordinary unions; callers and
            # existing replay tests may rely on them even when one image is
            # arithmetically redundant in another.
    return sp.simplify_logic(sp.Or(*vals))


def _scalar_value(result):
    from ..results import PropertyResult

    if isinstance(result, PropertyResult):
        return result.primary if result.primary is not None else None
    return result


def combine_piecewise_sign(values):
    vals = tuple(_scalar_value(v) for v in values)
    if not vals or any(v is None or v is Sign.UNKNOWN for v in vals):
        return Sign.UNKNOWN
    if all(v is Sign.ZERO for v in vals):
        return Sign.ZERO
    if all(v is Sign.POSITIVE for v in vals):
        return Sign.POSITIVE
    if all(v is Sign.NEGATIVE for v in vals):
        return Sign.NEGATIVE
    nonneg = {Sign.POSITIVE, Sign.NONNEGATIVE, Sign.ZERO}
    nonpos = {Sign.NEGATIVE, Sign.NONPOSITIVE, Sign.ZERO}
    if all(v in nonneg for v in vals):
        return Sign.NONNEGATIVE
    if all(v in nonpos for v in vals):
        return Sign.NONPOSITIVE
    if any(v in {Sign.POSITIVE, Sign.NONNEGATIVE} for v in vals) and any(
        v in {Sign.NEGATIVE, Sign.NONPOSITIVE} for v in vals
    ):
        return Sign.SIGN_CHANGING
    return Sign.UNKNOWN


def combine_piecewise_boundedness(values):
    vals = tuple(_scalar_value(v) for v in values)
    if not vals or any(v is None or v is Boundedness.UNKNOWN for v in vals):
        return Boundedness.UNKNOWN
    upper = {Boundedness.BOUNDED, Boundedness.BOUNDED_ABOVE}
    lower = {Boundedness.BOUNDED, Boundedness.BOUNDED_BELOW}
    above, below = all(v in upper for v in vals), all(v in lower for v in vals)
    if above and below:
        return Boundedness.BOUNDED
    if above:
        return Boundedness.BOUNDED_ABOVE
    if below:
        return Boundedness.BOUNDED_BELOW
    return Boundedness.UNBOUNDED


def combine_piecewise_monotonicity(values):
    vals = tuple(_scalar_value(v) for v in values)
    if not vals or any(v is None for v in vals):
        return Monotonicity.UNKNOWN
    inc = {Monotonicity.STRICTLY_INCREASING, Monotonicity.NONDECREASING}
    dec = {Monotonicity.STRICTLY_DECREASING, Monotonicity.NONINCREASING}
    if all(v is Monotonicity.CONSTANT for v in vals):
        return Monotonicity.CONSTANT
    if all(v in inc | {Monotonicity.CONSTANT} for v in vals):
        return Monotonicity.NONDECREASING
    if all(v in dec | {Monotonicity.CONSTANT} for v in vals):
        return Monotonicity.NONINCREASING
    if any(v in inc for v in vals) and any(v in dec for v in vals):
        return Monotonicity.NOT_MONOTONE
    return Monotonicity.UNKNOWN
