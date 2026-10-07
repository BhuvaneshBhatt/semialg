from __future__ import annotations

from builtins import range as builtin_range

import sympy as sp
from sympy import S

from .. import _semialg
from ..exactness import certify_identity_equal, certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import ConditionalExpression, Monotonicity, PropertyResult, TruthValue
from ..stratification import (
    combine_truth_results,
    exact_set_equal,
    guarded_result,
    normalize_cases,
)
from ..transcendental.periodic import positive_period
from .common import (
    combine_domain_condition,
    condition_to_set,
    is_semialgebraic_condition,
    is_semialgebraic_expression,
    normalize_variables,
)
from .domain import domain as _domain
from .monotonic_certificates import (
    injectivity_from_partition,
    monotonic_partition,
)
from .monotonicity import monotonicity
from .range import range


def _combined_domain(expr, var, domain):
    return combine_domain_condition(_domain(expr, var), domain, var)


def _semialg_mapping_result(expr, var, domain, *, codomain=True):
    """Return semialg mapping metadata when the whole problem is semialgebraic."""
    if not (
        is_semialgebraic_expression(expr, var)
        and is_semialgebraic_condition(domain, var)
    ):
        return None
    # semialg's mapping-properties query is currently nonparametric. Preserve
    # funcprops' guarded parameter strata
    # coefficients as part of an ordinary mapping query.
    problem_symbols = sp.sympify(expr).free_symbols | sp.sympify(domain).free_symbols
    if problem_symbols - {var}:
        return None
    y = sp.Dummy("_semialg_y", real=True)
    target = True
    if codomain not in (S.Reals, True, None, S.true):
        if isinstance(codomain, sp.Set):
            try:
                target = codomain.as_relational(y)
            except (TypeError, ValueError, NotImplementedError):
                return None
        else:
            target = sp.sympify(codomain).xreplace({var: y})
        if not is_semialgebraic_condition(target, y):
            return None
    try:
        return _semialg.function_mapping_properties(
            expr, var, domain=domain, codomain=target, image_variable=y
        )
    except (_semialg.SemialgUnavailable, *SYMBOLIC_QUERY_ERRORS):
        return None


def _truth_from_backend(value):
    return _semialg.normalize_truth(value)


def _domain_at_most_one_point(dom, var):
    ds = condition_to_set(dom, var)
    if ds is S.EmptySet:
        return True
    if isinstance(ds, sp.FiniteSet):
        return len(ds) <= 1
    return False


def _injectivity_from_monotonicity(m, dom, var):
    if m in (Monotonicity.STRICTLY_INCREASING, Monotonicity.STRICTLY_DECREASING):
        return TruthValue.TRUE
    if m is Monotonicity.CONSTANT:
        return (
            TruthValue.TRUE if _domain_at_most_one_point(dom, var) else TruthValue.FALSE
        )
    return TruthValue.UNKNOWN


def _multivariate_affine_map(expr, variables):
    """Return (A,b) for an affine map, else None.

    ``expr`` may be scalar or a tuple/list/Matrix. Coefficients may contain
    parameters but must not depend on the mapped variables.
    """
    vv = tuple(variables)
    if isinstance(expr, sp.MatrixBase):
        outs = tuple(expr)
    elif isinstance(expr, (tuple, list, sp.Tuple)):
        outs = tuple(map(sp.sympify, expr))
    else:
        outs = (sp.sympify(expr),)
    rows = []
    offsets = []
    for f in outs:
        try:
            row = [sp.simplify(sp.diff(f, v)) for v in vv]
            if any(a.has(*vv) for a in row):
                return None
            b = sp.simplify(f - sum(a * v for a, v in zip(row, vv, strict=False)))
            if b.has(*vv):
                return None
        except SYMBOLIC_QUERY_ERRORS:
            return None
        rows.append(row)
        offsets.append(b)
    return sp.Matrix(rows), sp.Matrix(offsets)


def _full_rank_condition(A, rank):
    """Return the exact symbolic condition that ``A`` has at least ``rank``."""
    from itertools import combinations

    m, n = A.shape
    if rank == 0:
        return S.true
    if rank > min(m, n):
        return S.false
    minors = []
    for rows in combinations(builtin_range(m), rank):
        for cols in combinations(builtin_range(n), rank):
            det = sp.simplify(A.extract(rows, cols).det())
            zero = certify_zero(det)
            if zero is False:
                return S.true
            if zero is None:
                minors.append(sp.Ne(det, 0))
    return sp.Or(*minors) if minors else S.false


def _rank_truth_result(A, required_rank, *, source):
    condition = sp.simplify(_full_rank_condition(A, required_rank))
    if condition is S.true:
        return TruthValue.TRUE
    if condition is S.false:
        return TruthValue.FALSE
    return normalize_cases(
        (
            ConditionalExpression(TruthValue.TRUE, condition),
            ConditionalExpression(TruthValue.FALSE, sp.Not(condition)),
        ),
        source=source,
    )


def _multivariate_mapping_truth(
    expr, variables, *, mode, domain=True, codomain=S.Reals
):
    # Exact support currently covers the natural full-real affine map. More
    # complicated restricted multivariate domains must not be guessed.
    if domain not in (True, None, S.true, S.Reals):
        return TruthValue.UNKNOWN
    aff = _multivariate_affine_map(expr, variables)
    if aff is None:
        return TruthValue.UNKNOWN
    A, _ = aff
    m, n = A.shape
    if mode == "injective":
        return _rank_truth_result(A, n, source="affine-rank-injectivity")
    if mode == "surjective":
        if codomain not in (S.Reals, True, None):
            return TruthValue.UNKNOWN
        return _rank_truth_result(A, m, source="affine-rank-surjectivity")
    if mode == "bijective":
        if m != n:
            return TruthValue.FALSE
        return _rank_truth_result(A, n, source="affine-rank-bijectivity")
    return TruthValue.UNKNOWN


def _affine_output_transform(expr, var):
    """Return ``(scale, inner, offset)`` for a nontrivial affine output map."""
    offset, dependent = expr.as_independent(var, as_Add=True)
    scale, inner = dependent.as_independent(var, as_Add=False)
    if inner == expr or (inner == dependent and scale == 1 and offset == 0):
        return None
    if scale.has(var) or offset.has(var) or not inner.has(var):
        return None
    return sp.simplify(scale), inner, sp.simplify(offset)


def _injective_structural_closure(expr, var, domain, partition):
    """Apply exact injectivity-preserving composition and affine-output rules."""
    transform = _affine_output_transform(expr, var)
    if transform is not None:
        scale, inner, _ = transform
        scale_zero = certify_zero(scale)
        if scale_zero is False:
            return injective(inner, var, domain=domain, partition=partition)
        if scale_zero is True:
            return (
                TruthValue.TRUE
                if _domain_at_most_one_point(domain, var)
                else TruthValue.FALSE
            )
        if scale.free_symbols:
            nonzero = injective(inner, var, domain=domain, partition=partition)
            zero = (
                TruthValue.TRUE
                if _domain_at_most_one_point(domain, var)
                else TruthValue.FALSE
            )
            nz = guarded_result(
                nonzero, sp.Ne(scale, 0), source="affine-output-injectivity"
            )
            z = guarded_result(
                zero, sp.Eq(scale, 0), source="affine-output-injectivity"
            )
            return normalize_cases(
                (*nz.cases, *z.cases), source="affine-output-injectivity"
            )
    # These outer real functions are strictly increasing on their natural real
    # domains, so composition preserves and reflects injectivity of the inner map.
    if (
        expr.func in {sp.exp, sp.log, sp.sinh, sp.asinh, sp.atan}
        and len(expr.args) == 1
    ):
        inner = expr.args[0]
        if inner != var and inner.has(var):
            return injective(inner, var, domain=domain, partition=partition)
    return TruthValue.UNKNOWN


def injective(expr, var, *, domain=True, partition=None):
    """Classify injectivity on the requested domain using exact witnesses.

    Scalar functions reuse monotonicity and monotonic partitions before trying
    exact equation solving. Multivariate cases use only supported structural
    criteria and otherwise return ``UNKNOWN``.
    """
    vv = normalize_variables(var)
    if len(vv) > 1:
        return _multivariate_mapping_truth(expr, vv, mode="injective", domain=domain)
    expr = sp.sympify(expr)
    var = vv[0]
    dom = _combined_domain(expr, var, domain)
    backend = _semialg_mapping_result(expr, var, dom)
    if backend is not None and backend.injective is not None:
        return _truth_from_backend(backend.injective)
    structural = _injective_structural_closure(expr, var, dom, partition)
    if structural is not TruthValue.UNKNOWN:
        return structural
    m = monotonicity(expr, var, domain=dom, partition=partition)
    if isinstance(m, PropertyResult):
        cases = []
        for c in m.cases:
            v = _injectivity_from_monotonicity(c.value, dom, var)
            if v is not TruthValue.UNKNOWN:
                cases.append(ConditionalExpression(v, c.condition))
        if cases:
            return normalize_cases(
                cases, source="monotonicity-injectivity-stratification"
            )
    else:
        direct = _injectivity_from_monotonicity(m, dom, var)
        if direct is not TruthValue.UNKNOWN:
            return direct

    # Reuse the exact monotonic partition
    # again.  Same-direction strict cells cover stationary inflections; a
    # certified direction reversal gives an exact noninjectivity proof.
    part = (
        partition
        if partition is not None
        and getattr(partition, "domain", None) == condition_to_set(dom, var)
        else monotonic_partition(expr, var, domain=dom, return_result=True)
    )
    pv = injectivity_from_partition(part)
    if pv is not TruthValue.UNKNOWN:
        return pv
    if not expr.has(var):
        return (
            TruthValue.TRUE if _domain_at_most_one_point(dom, var) else TruthValue.FALSE
        )

    # Periodicity is an exact disproof only when the requested interval contains
    # a pair x,x+p.  Merely having width == p is insufficient for a half-open
    # interval without additional continuity information.
    p = positive_period(expr, var)
    if p is not None:
        ds = condition_to_set(dom, var)
        if ds is S.Reals:
            return TruthValue.FALSE
        if isinstance(ds, sp.Interval) and ds.measure.is_number:
            d = sp.simplify(ds.measure - p)
            if d.is_positive is True:
                return TruthValue.FALSE
            if certify_zero(d) is True and not ds.left_open and not ds.right_open:
                return TruthValue.FALSE

    # Symmetry is a cheap exact disproof on symmetric domains, provided the
    # domain contains some nonzero point.
    ds = condition_to_set(dom, var)
    symmetric = ds is S.Reals or (
        isinstance(ds, sp.Interval)
        and certify_identity_equal(ds.start, -ds.end) is True
    )
    if symmetric:
        try:
            if certify_identity_equal(
                expr.xreplace({var: -var}), expr
            ) is True and expr.has(var):
                if ds is S.Reals or (isinstance(ds, sp.Interval) and ds.measure != 0):
                    return TruthValue.FALSE
        except SYMBOLIC_QUERY_ERRORS:
            pass
    return TruthValue.UNKNOWN


def _surjectivity_from_range(cond, y, target):
    rs = condition_to_set(cond, y)
    if rs is None:
        return TruthValue.UNKNOWN
    return exact_set_equal(rs, target)


def surjective(expr, var, *, domain=True, codomain=S.Reals):
    vv = normalize_variables(var)
    if len(vv) > 1:
        return _multivariate_mapping_truth(
            expr, vv, mode="surjective", domain=domain, codomain=codomain
        )
    expr = sp.sympify(expr)
    var = vv[0]
    dom = _combined_domain(expr, var, domain)
    backend = _semialg_mapping_result(expr, var, dom, codomain=codomain)
    if backend is not None and backend.surjective is not None:
        return _truth_from_backend(backend.surjective)
    y = sp.Dummy("_y", real=True)
    target = codomain if isinstance(codomain, sp.Set) else condition_to_set(codomain, y)
    if target is None:
        return TruthValue.UNKNOWN
    try:
        r = range(expr, var, y, constraints=dom)
    except SYMBOLIC_QUERY_ERRORS:
        return TruthValue.UNKNOWN
    if isinstance(r, PropertyResult):
        return normalize_cases(
            (
                ConditionalExpression(
                    _surjectivity_from_range(c.value, y, target), c.condition
                )
                for c in r.cases
            ),
            primary=_surjectivity_from_range(r.primary, y, target)
            if r.primary is not None
            else None,
            source="range-surjectivity-stratification",
        )
    return _surjectivity_from_range(r, y, target)


def bijective(expr, var, *, domain=True, codomain=S.Reals):
    vv = normalize_variables(var)
    if len(vv) > 1:
        return _multivariate_mapping_truth(
            expr, vv, mode="bijective", domain=domain, codomain=codomain
        )
    expr = sp.sympify(expr)
    var = vv[0]
    dom = _combined_domain(expr, var, domain)
    backend = _semialg_mapping_result(expr, var, dom, codomain=codomain)
    if backend is not None and backend.bijective is not None:
        return _truth_from_backend(backend.bijective)
    i = injective(expr, var, domain=dom)
    s = surjective(expr, var, domain=dom, codomain=codomain)
    return combine_truth_results(i, s, op="and", source="bijectivity-stratification")
