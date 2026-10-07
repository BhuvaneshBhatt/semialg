from __future__ import annotations

import sympy as sp
from sympy import S

from .. import _semialg
from ..engine import evaluate_property
from ..exactness import certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import ConditionalExpression, PropertyResult
from ..stratification import (
    evaluate_parameter_strata,
    normalize_cases,
    semantic_condition,
)
from ..transcendental.algebraize import algebraize
from .common import (
    affine_real_parameter_guard,
    combine_domain_condition,
    condition_to_set,
    is_full_natural_domain,
    set_to_condition,
)
from .domain import domain
from .monotonic_certificates import monotonic_partition, range_from_monotonic_partition
from .special_functions import (
    parameter_guard,
    positive_axis_range_strata,
    special_function_range_strata,
)


def _lambda_condition(value, y):
    """Instantiate a one-argument catalogue condition at ``y``."""
    if isinstance(value, sp.Lambda) and len(value.variables) == 1:
        return sp.sympify(value.expr).xreplace({value.variables[0]: y})
    return None


def _refine_interval_attainment(expr, var, cond, y):
    s = condition_to_set(cond, y)
    if not isinstance(s, sp.Interval):
        return cond
    natural = domain(expr, var)
    ds = condition_to_set(natural, var)
    if ds is None:
        return cond
    left_open = s.left_open
    right_open = s.right_open
    if s.start not in (S.NegativeInfinity, S.Infinity):
        try:
            sol = sp.solveset(sp.Eq(expr, s.start), var, domain=ds)
            if sol is S.EmptySet:
                left_open = True
            elif isinstance(sol, (sp.FiniteSet, sp.ImageSet, sp.Union)):
                left_open = False
        except SYMBOLIC_QUERY_ERRORS:
            pass
    if s.end not in (S.NegativeInfinity, S.Infinity):
        try:
            sol = sp.solveset(sp.Eq(expr, s.end), var, domain=ds)
            if sol is S.EmptySet:
                right_open = True
            elif isinstance(sol, (sp.FiniteSet, sp.ImageSet, sp.Union)):
                right_open = False
        except SYMBOLIC_QUERY_ERRORS:
            pass
    return set_to_condition(
        sp.Interval(s.start, s.end, left_open=left_open, right_open=right_open), y
    )


def _registry_range(expr, var, y):
    r = evaluate_property(expr, "Range", var=var, mode="real")
    pguard = parameter_guard(expr, var, "Range")
    if r.primary is not None:
        cond = _lambda_condition(r.primary, y)
        if cond is None:
            return None
        cond = _refine_interval_attainment(expr, var, cond, y)
        if pguard is S.true:
            return cond
        return normalize_cases(
            (ConditionalExpression(cond, pguard),),
            source="registry-range-stratification",
        )
    if r.cases:
        cases = []
        for c in r.cases:
            cond = _lambda_condition(c.value, y)
            if cond is not None:
                cond = _refine_interval_attainment(expr, var, cond, y)
                cases.append(
                    ConditionalExpression(
                        cond, sp.And(semantic_condition(c.condition), pguard)
                    )
                )
        if cases:
            return normalize_cases(cases, source="registry-range-stratification")
    return None


def _monotone_image_interval(outer, X, inner_set):
    from ..results import Monotonicity
    from .monotonicity import monotonicity

    dom = set_to_condition(inner_set, X)
    m = monotonicity(outer, X, domain=dom)
    inc = m in (Monotonicity.STRICTLY_INCREASING, Monotonicity.NONDECREASING)
    dec = m in (Monotonicity.STRICTLY_DECREASING, Monotonicity.NONINCREASING)
    if not (inc or dec):
        return None

    def image_endpoint(v, is_left, is_open):
        if v is S.NegativeInfinity:
            return sp.limit(outer, X, -sp.oo)
        if v is S.Infinity:
            return sp.limit(outer, X, sp.oo)
        if is_open:
            try:
                return sp.limit(outer, X, v, dir="+" if is_left else "-")
            except SYMBOLIC_QUERY_ERRORS:
                pass
        return sp.simplify(outer.subs(X, v))

    a = image_endpoint(inner_set.start, True, inner_set.left_open)
    b = image_endpoint(inner_set.end, False, inner_set.right_open)
    if inc:
        lo, hi = a, b
        lo_open = inner_set.left_open or inner_set.start in (
            S.NegativeInfinity,
            S.Infinity,
        )
        hi_open = inner_set.right_open or inner_set.end in (
            S.NegativeInfinity,
            S.Infinity,
        )
    else:
        lo, hi = b, a
        lo_open = inner_set.right_open or inner_set.end in (
            S.NegativeInfinity,
            S.Infinity,
        )
        hi_open = inner_set.left_open or inner_set.start in (
            S.NegativeInfinity,
            S.Infinity,
        )
    if lo is S.NegativeInfinity:
        lo_open = True
    if hi is S.Infinity:
        hi_open = True
    return sp.Interval(lo, hi, left_open=lo_open, right_open=hi_open)


def _composition_range(expr, var, y, domain_condition):
    varying = [(i, a) for i, a in enumerate(getattr(expr, "args", ())) if a.has(var)]
    if len(varying) != 1 or expr.func in (sp.Add, sp.Mul, sp.Pow, sp.Piecewise):
        return None
    i, inner = varying[0]
    # A primitive f(x) is not a nontrivial composition.  Recursing through
    # range(x, ...) and then asking monotonicity of the same primitive
    # can manufacture bogus empty/False images for oscillatory functions.
    if inner == var:
        return None
    u = sp.Dummy("u", real=True)
    rc = range(inner, var, u, constraints=domain_condition)
    if isinstance(rc, PropertyResult):
        return None
    inner_set = condition_to_set(rc, u)
    if not isinstance(inner_set, sp.Interval):
        return None
    X = sp.Dummy("X", real=True)
    args = list(expr.args)
    args[i] = X
    outer = expr.func(*args)
    out = _monotone_image_interval(outer, X, inner_set)
    return set_to_condition(out, y) if out is not None else None


def _scale_interval(s, A, B, sign):
    def finite_image(v):
        if v is S.NegativeInfinity:
            return S.NegativeInfinity if sign > 0 else S.Infinity
        if v is S.Infinity:
            return S.Infinity if sign > 0 else S.NegativeInfinity
        return sp.simplify(A * v + B)

    if sign > 0:
        return sp.Interval(
            finite_image(s.start),
            finite_image(s.end),
            left_open=s.left_open,
            right_open=s.right_open,
        )
    return sp.Interval(
        finite_image(s.end),
        finite_image(s.start),
        left_open=s.right_open,
        right_open=s.left_open,
    )


def _affine_transform(expr, var, y, domain):
    """Transform an exact inner range through ``A*g(x) + B``.

    Parameter-dependent coefficient signs are preserved as conditional strata.
    """
    # Identify A*g(x)+B where A and B are independent of the input variable.
    terms = sp.Add.make_args(sp.expand(expr))
    const = sum((t for t in terms if not t.has(var)), S.Zero)
    dep = sp.simplify(expr - const)
    coeff, rest = dep.as_independent(var, as_Add=False)
    if coeff.has(var) or not rest.has(var) or rest == expr:
        return None
    u = sp.Dummy("u", real=True)
    inner = range(rest, var, u, constraints=domain)

    def transform(cond, extra=S.true):
        ss = condition_to_set(cond, u)
        if not isinstance(ss, sp.Interval):
            return None
        A = sp.sympify(coeff)
        B = sp.sympify(const)
        if A.is_positive:
            return set_to_condition(_scale_interval(ss, A, B, 1), y)
        if A.is_negative:
            return set_to_condition(_scale_interval(ss, A, B, -1), y)
        if certify_zero(A) is True:
            return sp.Eq(y, B)
        return normalize_cases(
            (
                ConditionalExpression(
                    set_to_condition(_scale_interval(ss, A, B, 1), y),
                    sp.And(extra, A > 0),
                ),
                ConditionalExpression(sp.Eq(y, B), sp.And(extra, sp.Eq(A, 0))),
                ConditionalExpression(
                    set_to_condition(_scale_interval(ss, A, B, -1), y),
                    sp.And(extra, A < 0),
                ),
            ),
            source="affine-range-stratification",
        )

    if isinstance(inner, PropertyResult):
        cases = []
        for c in inner.cases:
            tr = transform(c.value, c.condition)
            if isinstance(tr, PropertyResult):
                cases.extend(tr.cases)
            elif tr is not None:
                cases.append(ConditionalExpression(tr, c.condition))
        return (
            normalize_cases(cases, source="affine-range-stratification")
            if cases
            else None
        )
    return transform(inner)


def _parameterized_elementary_range(expr, var, y, domain_condition):
    """High-value exact parameter strata not handled by generic calculus."""
    ds = condition_to_set(domain_condition, var)
    try:
        positive_axis = (
            isinstance(ds, sp.Interval)
            and ds.start == 0
            and ds.left_open
            and ds.end is S.Infinity
        ) or sp.simplify_logic(
            sp.Implies(semantic_condition(domain_condition), var > 0)
        ) is S.true
    except SYMBOLIC_QUERY_ERRORS:
        positive_axis = False
    # x**a on x>0: every nonzero real exponent maps onto (0,oo), while
    # a=0 is the singleton 1.  This must precede global catalogue ranges.
    if positive_axis and expr.func is sp.Pow and expr.base == var:
        A = expr.exp
        if A.free_symbols:
            real_guard = S.true if A.is_real is True else sp.Contains(A, S.Reals)
            return normalize_cases(
                (
                    ConditionalExpression(y > 0, sp.And(real_guard, sp.Ne(A, 0))),
                    ConditionalExpression(sp.Eq(y, 1), sp.Eq(A, 0)),
                ),
                source="positive-power-range-stratification",
            )
    if ds is not S.Reals:
        return None
    # exp(a*x+b): nonzero real slope maps R bijectively to (0,oo); zero
    # slope is the singleton exp(b).
    if expr.func is sp.exp:
        affine = affine_real_parameter_guard(expr.args[0], var)
        if affine is not None:
            A, B, real_guard = affine
            if A.free_symbols or B.free_symbols or real_guard is not S.true:
                real_b = S.true if B.is_real is True else sp.Q.real(B)
                return normalize_cases(
                    (
                        ConditionalExpression(y > 0, sp.And(real_guard, sp.Ne(A, 0))),
                        ConditionalExpression(
                            sp.Eq(y, sp.exp(B)), sp.And(real_b, sp.Eq(A, 0))
                        ),
                    ),
                    source="exp-affine-range-stratification",
                )
    return None


def _sympy_range(expr, var, y, domain_condition):
    ds = condition_to_set(domain_condition, var)
    if ds is None:
        return None
    try:
        rs = sp.calculus.util.function_range(expr, var, ds)
    except SYMBOLIC_QUERY_ERRORS:
        return None
    if isinstance(rs, sp.Set) and not rs.has(var):
        return set_to_condition(rs, y)
    return None


def _early_range_stage(expr, var, y, domain_condition, constraints):
    """Run cheap parameter, catalogue, and positive-axis range rules."""
    pr = _parameterized_elementary_range(expr, var, y, domain_condition)
    if pr is not None:
        return pr
    sf_domain = domain_condition if constraints is True else sp.sympify(constraints)
    sfacts = special_function_range_strata(expr, var, y, sf_domain)
    if sfacts:
        out = normalize_cases(
            (ConditionalExpression(v, g) for v, g, *_ in sfacts),
            source="dlmf-special-function-range-stratification",
            evidence=tuple(src for *_, src in sfacts),
        )
        return out.primary if out.primary is not None else out
    ds0 = condition_to_set(domain_condition, var)
    try:
        positive_axis = (
            isinstance(ds0, sp.Interval)
            and ds0.start == 0
            and ds0.left_open
            and ds0.end is S.Infinity
        ) or sp.simplify_logic(
            sp.Implies(semantic_condition(domain_condition), var > 0)
        ) is S.true
    except SYMBOLIC_QUERY_ERRORS:
        positive_axis = False
    if positive_axis:
        facts = positive_axis_range_strata(expr, var, y)
        if facts:
            out = normalize_cases(
                (ConditionalExpression(v, g) for v, g, *_ in facts),
                source="special-function-range-stratification",
                evidence=tuple(src for *_, src in facts),
            )
            return out.primary if out.primary is not None else out
    return None


def _semialgebraic_range_stage(expr, var, y, domain_condition, method, kwargs):
    """Run structural then exact semialgebraic range analysis."""
    if not _semialg.eligible_problem(expr, var, domain=domain_condition):
        return None
    structural = _affine_transform(expr, var, y, domain_condition)
    if structural is not None:
        return (
            structural
            if isinstance(structural, PropertyResult)
            else sp.simplify(structural)
        )
    try:
        return _semialg.range(
            expr,
            var,
            y,
            constraints=domain_condition,
            method="qe" if method == "auto" else method,
            **kwargs,
        )
    except SYMBOLIC_QUERY_ERRORS:
        return None


def _piecewise_range(expr, var, y, domain_condition, method, kwargs):
    from .piecewise import combine_piecewise_conditions, piecewise_decomposition

    decomp = piecewise_decomposition(expr, var, domain_condition)
    if not decomp.complete:
        return None
    images = []
    for branch in decomp.branches:
        image = _function_range_impl(
            branch.expression,
            var,
            y,
            constraints=branch.condition,
            method=method,
            partition=None,
            _domain_ready=True,
            **kwargs,
        )
        if isinstance(image, PropertyResult):
            return None
        images.append(image)
    return combine_piecewise_conditions(images, y)


def _function_range_impl(
    expr,
    var,
    y,
    *,
    constraints=True,
    method="auto",
    partition=None,
    _domain_ready=False,
    **kwargs,
):
    """Return the real range as a Boolean condition in ``y``.

    Strategy: primitive lookup, structural transformations, SymPy's univariate
    calculus algorithms, and semialg for the semialgebraic subset.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    y = sp.sympify(y)
    natural = domain(expr, var)
    domain_condition = (
        sp.sympify(constraints)
        if _domain_ready
        else combine_domain_condition(natural, constraints, var)
    )
    if not _domain_ready:
        split = evaluate_parameter_strata(
            domain_condition,
            var,
            lambda branch_domain: _function_range_impl(
                expr,
                var,
                y,
                constraints=branch_domain,
                method=method,
                partition=None,
                _domain_ready=True,
                **kwargs,
            ),
            source="range-parameter-stratification",
        )
        if split is not None:
            return split
    if not expr.has(var):
        return sp.Eq(y, expr)

    if isinstance(expr, sp.Piecewise):
        merged = _piecewise_range(expr, var, y, domain_condition, method, kwargs)
        if merged is not None:
            return merged

    early = _early_range_stage(expr, var, y, domain_condition, constraints)
    if early is not None:
        return early

    # Preserve parameter strata before delegating semialgebraic elimination.
    # A backend QE call may quantify/treat free coefficients in a way that
    # collapses guarded image cases (for example a*x on x>0) to a tautology.
    backend = _semialgebraic_range_stage(expr, var, y, domain_condition, method, kwargs)
    if backend is not None:
        return backend

    # Exact transcendental algebraization can turn the range problem into a
    # semialgebraic one without changing the output variable. For example,
    # exp(2*x)-3*exp(x)+2 becomes u**2-3*u+2 on u>0.
    ar = algebraize(expr, var)
    if ar.changed and ar.exact and len(ar.substitutions) == 1:
        u, _ = ar.substitutions[0]
        # Preserve openness of the auxiliary image (notably u=exp(x)>0).
        # SymPy's univariate calculus is tried first because it represents
        # endpoint attainment explicitly; semialg remains the exact fallback.
        try:
            us = condition_to_set(ar.constraints, u)
            if us is not None:
                rs = sp.calculus.util.function_range(ar.expression, u, us)
                if isinstance(rs, sp.Set):
                    return set_to_condition(rs, y)
        except SYMBOLIC_QUERY_ERRORS:
            pass
        try:
            raw = _semialg.range(
                ar.expression,
                u,
                y,
                constraints=ar.constraints,
                method="qe" if method == "auto" else method,
                **kwargs,
            )
            # Refine finite interval endpoints against the actual auxiliary
            # domain so a backend closure approximation cannot turn > into >=.
            rs = condition_to_set(raw, y)
            us = condition_to_set(ar.constraints, u)
            if isinstance(rs, sp.Interval) and us is not None:
                lo_open, hi_open = rs.left_open, rs.right_open
                if rs.start not in (S.NegativeInfinity, S.Infinity):
                    sol = sp.solveset(sp.Eq(ar.expression, rs.start), u, domain=us)
                    if sol is S.EmptySet:
                        lo_open = True
                if rs.end not in (S.NegativeInfinity, S.Infinity):
                    sol = sp.solveset(sp.Eq(ar.expression, rs.end), u, domain=us)
                    if sol is S.EmptySet:
                        hi_open = True
                return set_to_condition(
                    sp.Interval(
                        rs.start, rs.end, left_open=lo_open, right_open=hi_open
                    ),
                    y,
                )
            return raw
        except SYMBOLIC_QUERY_ERRORS:
            pass

    # Monotone composition can turn a known inner interval into the exact
    # outer range (e.g. exp(sin(x)) -> [exp(-1), exp(1)]).
    r = _composition_range(expr, var, y, domain_condition)
    if r is not None:
        return r if isinstance(r, PropertyResult) else sp.simplify(r)

    # Primitive catalogue ranges are exact and cheapest.
    if constraints is True or is_full_natural_domain(domain_condition, natural):
        r = _registry_range(expr, var, y)
        if r is not None:
            return r if isinstance(r, PropertyResult) else sp.simplify(r)

    ds = condition_to_set(domain_condition, var)
    part = (
        partition
        if partition is not None and getattr(partition, "domain", None) == ds
        else None
    )
    # Indexed special-function geometry is cheap to recognize and expensive to
    # rediscover through generic calculus, so build the shared partition first
    # for those families.  For generic expressions keep SymPy's successful fast
    # path ahead of a potentially redundant derivative-root solve.
    if isinstance(ds, sp.Interval):
        try:
            from ..transcendental.zero_geometry import zero_geometry

            if part is None and zero_geometry(expr, var, derivative_order=1, domain=ds):
                part = monotonic_partition(expr, var, domain=ds, return_result=True)
        except SYMBOLIC_QUERY_ERRORS:
            pass
    if part is not None:
        prange = range_from_monotonic_partition(expr, var, part)
        if prange.complete and isinstance(prange.image, sp.Set):
            return set_to_condition(prange.image, y)

    r = _sympy_range(expr, var, y, domain_condition)
    if r is not None:
        return sp.simplify(r)

    if part is None and isinstance(ds, sp.Interval):
        part = monotonic_partition(expr, var, domain=ds, return_result=True)
        prange = range_from_monotonic_partition(expr, var, part)
        if prange.complete and isinstance(prange.image, sp.Set):
            return set_to_condition(prange.image, y)

    # Indexed extrema consume the same partition certificate, avoiding a second
    # derivative-zero geometry pass.
    try:
        from ..transcendental.extrema import symbolic_extremum_range

        er = symbolic_extremum_range(
            expr,
            var,
            ds if isinstance(ds, sp.Interval) else domain_condition,
            partition=part,
        )
        if er is not None and er.complete:
            cond = er.as_condition(y)
            return cond if cond is not None else er
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return sp.Function("UnknownRange")(expr, var, y)


def range(expr, var, y, *, constraints=True, method="auto", partition=None, **kwargs):
    """Return the real range as a Boolean condition in ``y``."""
    return _function_range_impl(
        expr,
        var,
        y,
        constraints=constraints,
        method=method,
        partition=partition,
        _domain_ready=False,
        **kwargs,
    )
