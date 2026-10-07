"""Joint continuity, differentiability, smoothness, and analyticity analysis."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from .. import _semialg
from ..exactness import certify_identity_equal
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import Analyticity, Continuity, TruthValue
from ..singularities import (
    analytic,
    continuous,
    discontinuities,
)
from ..stratification import evaluate_parameter_strata
from .common import (
    combine_domain_condition,
    condition_to_set,
)
from .domain import domain as _domain


@dataclass(frozen=True)
class SmoothnessResult:
    """Certified regularity information on a requested real domain."""

    continuity: Continuity
    differentiability_order: int | sp.Expr | None
    smooth: bool | None
    analytic: Analyticity
    exceptional_loci: tuple = ()
    evidence: tuple = ()


def _exception_loci(records):
    """Return unique singularity conditions while preserving source order."""
    loci = []
    for record in records:
        condition = getattr(record, "condition", record)
        if condition not in loci:
            loci.append(condition)
    return tuple(loci)


def _locus_intersects_domain(locus, domain, var):
    """Conservatively decide whether a locus can meet a requested domain."""
    domain_set = condition_to_set(domain, var)
    locus_set = condition_to_set(locus, var)
    if domain_set is not None and locus_set is not None:
        try:
            overlap = sp.Intersection(domain_set, locus_set)
        except (TypeError, ValueError, NotImplementedError):
            overlap = None
        if overlap is not None:
            if overlap is S.EmptySet or overlap == S.EmptySet:
                return False
            is_empty = getattr(overlap, "is_empty", None)
            if is_empty is True:
                return False
            if is_empty is False:
                return True

    try:
        satisfiable = sp.satisfiable(sp.And(sp.sympify(domain), sp.sympify(locus)))
    except (TypeError, ValueError, NotImplementedError):
        return True
    return satisfiable is not False


def _domain_loci(records, domain, var):
    return tuple(
        locus
        for locus in _exception_loci(records)
        if _locus_intersects_domain(locus, domain, var)
    )


def _structurally_smooth(expr, var):
    """Recognize globally smooth expression trees without derivative probing.

    The rule set intentionally contains only closure operations that cannot
    introduce a real singularity: polynomial arithmetic, nonnegative integer
    powers, and composition with entire elementary functions.
    """
    expr = sp.sympify(expr)
    if not expr.has(var):
        return True
    if expr == var:
        return True
    if expr.is_Add or expr.is_Mul:
        return all(_structurally_smooth(arg, var) for arg in expr.args)
    if expr.is_Pow:
        base, exponent = expr.args
        return bool(exponent.is_Integer and exponent.is_nonnegative) and (
            _structurally_smooth(base, var)
        )
    entire_heads = {sp.exp, sp.sin, sp.cos, sp.sinh, sp.cosh, sp.erf}
    if expr.func in entire_heads and len(expr.args) == 1:
        return _structurally_smooth(expr.args[0], var)
    return False


def _smoothness_value(result):
    """Return the certified differentiability order.

    ``0`` denotes C^0 regularity, positive integers denote the highest
    certified finite C^k order, ``oo`` denotes C-infinity, and ``None``
    means that no C^k regularity classification was certified.
    """
    return result.differentiability_order


def smooth(expr, var, *, domain=True, max_order=6):
    """Return whether ``expr`` is certified C-infinity on ``domain``.

    This is the predicate companion to :func:`smoothness`, whose
    mathematical default is the certified differentiability order.
    """
    result = _function_smoothness_impl(
        expr, var, domain=domain, max_order=max_order, return_result=True
    )
    if isinstance(result, SmoothnessResult):
        if result.smooth is True:
            return TruthValue.TRUE
        if result.smooth is False:
            return TruthValue.FALSE
    return TruthValue.UNKNOWN


def _smoothness_output(result, *, return_result):
    return result if return_result else _smoothness_value(result)


def _piecewise_smoothness(expr, var, requested_domain, max_order):
    """Certify finite-order seam regularity without differentiating Piecewise itself."""
    from .piecewise import certify_piecewise_boundaries, piecewise_decomposition

    decomp = piecewise_decomposition(expr, var, requested_domain)
    if not decomp.complete or not decomp.branches:
        return None

    branch_results = []
    for branch in decomp.branches:
        result = _function_smoothness_impl(
            branch.expression,
            var,
            domain=branch.condition,
            max_order=max_order,
            return_result=True,
            _domain_ready=True,
        )
        if not isinstance(result, SmoothnessResult):
            return None
        branch_results.append(result)

    seams = certify_piecewise_boundaries(decomp, max_order=max_order)
    if any(c.continuous is False for c in seams):
        return SmoothnessResult(
            Continuity.DISCONTINUOUS,
            0,
            False,
            Analyticity.NOT_ANALYTIC,
            tuple(c.boundary.point for c in seams if c.continuous is False),
            (decomp, *seams),
        )
    if any(c.continuous is None for c in seams):
        continuity = Continuity.UNKNOWN
    else:
        continuity = (
            Continuity.CONTINUOUS
            if all(r.continuity is Continuity.CONTINUOUS for r in branch_results)
            else Continuity.UNKNOWN
        )

    finite_orders = []
    for r in branch_results:
        if r.differentiability_order is sp.oo or r.differentiability_order == sp.oo:
            finite_orders.append(max_order)
        elif isinstance(r.differentiability_order, int):
            finite_orders.append(r.differentiability_order)
        else:
            finite_orders.append(0)
    for seam in seams:
        finite_orders.append(seam.differentiability_order or 0)
    verified = min(finite_orders) if finite_orders else 0

    # A witnessed seam mismatch proves non-smoothness.  Matching finitely many
    # derivatives never proves C-infinity.  The one safe nontrivial exception is
    # when adjacent branch formulas are symbolically identical, in which case
    # the Piecewise syntax is artificial and the common expression can be used.
    mismatch = any(
        seam.continuous is False
        or (
            seam.continuous is True
            and seam.differentiability_order is not None
            and seam.differentiability_order < max_order
        )
        for seam in seams
    )
    same_formula = all(
        certify_identity_equal(branch.expression, decomp.branches[0].expression) is True
        for branch in decomp.branches[1:]
    )
    if same_formula:
        return _function_smoothness_impl(
            decomp.branches[0].expression,
            var,
            domain=requested_domain,
            max_order=max_order,
            return_result=True,
            _domain_ready=True,
        )
    analytic = (
        Analyticity.NOT_ANALYTIC
        if seams
        else (
            Analyticity.ANALYTIC
            if all(r.analytic is Analyticity.ANALYTIC for r in branch_results)
            else Analyticity.UNKNOWN
        )
    )
    return SmoothnessResult(
        continuity,
        verified,
        False if mismatch else None,
        analytic,
        tuple(
            c.boundary.point
            for c in seams
            if c.continuous is not True or (c.differentiability_order or 0) < max_order
        ),
        (decomp, *branch_results, *seams),
    )


def _function_smoothness_impl(
    expr,
    var,
    *,
    domain=True,
    max_order=6,
    return_result=False,
    _domain_ready=False,
):
    """Analyze regularity up to ``max_order`` on the requested real domain.

    Finite derivative probing never certifies :math:`C^\\infty`.  Smoothness is
    true only when analyticity proves it; it is false only after a witnessed
    finite-order failure, and otherwise remains unknown.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if isinstance(max_order, bool) or not isinstance(max_order, int) or max_order < 1:
        raise ValueError("max_order must be a positive integer")

    natural_domain = _domain(expr, var)
    requested_domain = (
        sp.sympify(domain)
        if _domain_ready
        else combine_domain_condition(natural_domain, domain, var)
    )
    if not _domain_ready:
        split = evaluate_parameter_strata(
            requested_domain,
            var,
            lambda branch_domain: _function_smoothness_impl(
                expr,
                var,
                domain=branch_domain,
                max_order=max_order,
                return_result=return_result,
                _domain_ready=True,
            ),
            source="smoothness-parameter-stratification",
        )
        if split is not None:
            return split
    if isinstance(expr, sp.Piecewise):
        piecewise = _piecewise_smoothness(expr, var, requested_domain, max_order)
        if piecewise is not None:
            return _smoothness_output(piecewise, return_result=return_result)
    analyticity = analytic(expr, var, domain=requested_domain)
    if analyticity is Analyticity.ANALYTIC:
        continuity = continuous(expr, var, domain=requested_domain)
        exceptions = list(
            _domain_loci(
                discontinuities(expr, var, return_result=True),
                requested_domain,
                var,
            )
        )
        return _smoothness_output(
            SmoothnessResult(
                continuity,
                sp.oo,
                True,
                analyticity,
                tuple(exceptions),
                ("analytic-implies-smooth",),
            ),
            return_result=return_result,
        )

    backend = None
    if _semialg.eligible_problem(expr, var, domain=requested_domain):
        try:
            backend = _semialg.smoothness(
                expr, var, domain=requested_domain, max_order=max_order
            )
            # Negative/finite regularity results are safe to use immediately.
            # A backend claim of C-infinity is validated below through the
            # requested finite derivative order before it can strengthen the
            # public result.
            if backend.smooth is not True:
                continuity = (
                    Continuity.CONTINUOUS
                    if backend.continuous is True
                    else Continuity.DISCONTINUOUS
                    if backend.continuous is False
                    else Continuity.UNKNOWN
                )
                exceptions = [
                    backend.continuity_exceptions,
                    *backend.derivative_exceptions.values(),
                ]
                exceptions = [
                    locus for locus in exceptions if locus not in (S.false, False)
                ]
                return _smoothness_output(
                    SmoothnessResult(
                        continuity,
                        backend.differentiability_order,
                        backend.smooth,
                        analyticity,
                        tuple(exceptions),
                        (("semialg", backend),),
                    ),
                    return_result=return_result,
                )
        except (_semialg.SemialgUnavailable, *SYMBOLIC_QUERY_ERRORS):
            backend = None

    continuity = continuous(expr, var, domain=requested_domain)
    exceptions = list(
        _domain_loci(
            discontinuities(expr, var, return_result=True),
            requested_domain,
            var,
        )
    )

    if continuity is Continuity.CONTINUOUS and _structurally_smooth(expr, var):
        return _smoothness_output(
            SmoothnessResult(
                continuity,
                sp.oo,
                True,
                analyticity,
                tuple(exceptions),
                ("smooth-structural-closure",),
            ),
            return_result=return_result,
        )

    order = 0 if continuity is Continuity.CONTINUOUS else None
    failure_proved = False
    domain_set = condition_to_set(requested_domain, var)
    if domain_set is None:
        return _smoothness_output(
            SmoothnessResult(
                continuity, order, None, analyticity, tuple(exceptions), ()
            ),
            return_result=return_result,
        )

    from ..proof_cache import property_proof_context

    context = property_proof_context(expr, var, domain=requested_domain)
    for derivative_order in range(1, max_order + 1):
        try:
            derivative = context.derivative(derivative_order)
        except (TypeError, ValueError, NotImplementedError):
            break

        loci = _domain_loci(
            discontinuities(derivative, var, return_result=True),
            requested_domain,
            var,
        )
        for locus in loci:
            if locus not in exceptions:
                exceptions.append(locus)
        if loci:
            failure_proved = True
            break

        try:
            continuous_domain = sp.calculus.util.continuous_domain(
                derivative, var, domain_set
            )
        except (TypeError, ValueError, NotImplementedError):
            # An unsupported continuity query is uncertainty, not evidence of
            # either continuity or discontinuity.
            break

        if continuous_domain == domain_set:
            order = derivative_order
            continue

        missing = sp.Complement(domain_set, continuous_domain)
        if missing is not S.EmptySet and missing != S.EmptySet:
            failure_proved = True
        break

    smooth = False if failure_proved else None
    evidence = (("verified-derivatives-through", order),) if order is not None else ()
    if (
        not failure_proved
        and backend is not None
        and backend.smooth is True
        and order == max_order
    ):
        order = sp.oo
        smooth = True
        evidence = (*evidence, ("semialg", backend))
    return _smoothness_output(
        SmoothnessResult(
            continuity,
            order,
            smooth,
            analyticity,
            tuple(exceptions),
            evidence,
        ),
        return_result=return_result,
    )


def smoothness(expr, var, *, domain=True, max_order=6, return_result=False):
    """Return the certified C^k order, or detailed regularity information.

    The mathematical default is ``0`` for C^0, a positive integer for a
    certified finite C^k order, ``oo`` for C-infinity, or ``None`` when no
    regularity order is certified.  Use ``return_result=True`` for the full
    :class:`SmoothnessResult`.
    """
    return _function_smoothness_impl(
        expr,
        var,
        domain=domain,
        max_order=max_order,
        return_result=return_result,
        _domain_ready=False,
    )
