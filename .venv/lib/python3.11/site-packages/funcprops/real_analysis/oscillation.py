from __future__ import annotations

import contextlib
from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import Oscillation
from ..stratification import evaluate_parameter_strata
from ..transcendental.periodic import positive_period
from .common import combine_domain_condition, condition_to_set
from .domain import domain as _domain
from .monotonic_certificates import MonotonicPartition, monotonic_partition


@dataclass(frozen=True)
class OscillationCertificate:
    classification: Oscillation
    domain: sp.Set
    partition: MonotonicPartition | None = None
    turning_count: object = None
    source: str = ""
    evidence: tuple = ()


def _function_oscillation_impl(
    expr,
    var,
    *,
    domain=True,
    return_result=False,
    partition=None,
    _domain_ready=False,
):
    """Classify real oscillation from the shared monotonic partition.

    ``OSCILLATORY`` means a certified infinite sequence of monotonicity
    reversals (or a nonconstant exact periodic function on an unbounded
    interval). ``FINITELY_OSCILLATORY`` means a complete finite partition with
    at least one direction reversal.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    natural = _domain(expr, var)
    cond = (
        sp.sympify(domain)
        if _domain_ready
        else combine_domain_condition(natural, domain, var)
    )
    with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
        cond = sp.simplify(cond)
    if not _domain_ready:
        split = evaluate_parameter_strata(
            cond,
            var,
            lambda branch_domain: _function_oscillation_impl(
                expr,
                var,
                domain=branch_domain,
                return_result=return_result,
                partition=None,
                _domain_ready=True,
            ),
            source="oscillation-parameter-stratification",
        )
        if split is not None:
            return split
    ds = condition_to_set(cond, var)
    if not isinstance(ds, sp.Set):
        out = OscillationCertificate(
            Oscillation.UNKNOWN, S.Reals, None, None, "unknown-domain"
        )
        return out if return_result else out.classification
    if not expr.has(var):
        out = OscillationCertificate(
            Oscillation.NONOSCILLATORY, ds, None, 0, "constant"
        )
        return out if return_result else out.classification
    part = (
        partition
        if partition is not None and getattr(partition, "domain", None) == ds
        else monotonic_partition(expr, var, domain=ds, return_result=True)
    )
    if part.indexed is not None and part.indexed.boundary_geometry:
        try:
            from ..transcendental.zero_geometry import geometry_support_inside_domain

            infinite_inside = any(
                g.count is sp.oo and geometry_support_inside_domain(g, ds)
                for g in part.indexed.boundary_geometry
            )
        except SYMBOLIC_QUERY_ERRORS:
            infinite_inside = False
        if infinite_inside:
            out = OscillationCertificate(
                Oscillation.OSCILLATORY,
                ds,
                part,
                sp.oo,
                "indexed-monotonic-partition",
                (part,),
            )
            return out if return_result else out.classification
    if part.complete and part.cells:
        turns = part.direction_changes
        cls = Oscillation.FINITELY_OSCILLATORY if turns else Oscillation.NONOSCILLATORY
        out = OscillationCertificate(
            cls, ds, part, turns, "finite-monotonic-partition", (part,)
        )
        return out if return_result else out.classification
    # A certified global monotonic classification is a cheap exact proof of
    # nonoscillation when critical-point partitioning itself is incomplete.
    try:
        from ..results import Monotonicity
        from .monotonicity import monotonicity

        m = monotonicity(expr, var, domain=cond, partition=part)
        if m in (
            Monotonicity.STRICTLY_INCREASING,
            Monotonicity.NONDECREASING,
            Monotonicity.CONSTANT,
            Monotonicity.NONINCREASING,
            Monotonicity.STRICTLY_DECREASING,
        ):
            out = OscillationCertificate(
                Oscillation.NONOSCILLATORY,
                ds,
                part,
                0,
                "global-monotonicity",
                (m, part),
            )
            return out if return_result else out.classification
    except SYMBOLIC_QUERY_ERRORS:
        pass
    p = positive_period(expr, var)
    if (
        p is not None
        and isinstance(ds, sp.Interval)
        and (ds.start is -sp.oo or ds.end is sp.oo)
    ):
        out = OscillationCertificate(
            Oscillation.OSCILLATORY, ds, part, sp.oo, "exact-periodicity", (p, part)
        )
        return out if return_result else out.classification
    out = OscillationCertificate(
        Oscillation.UNKNOWN, ds, part, None, "incomplete-monotonic-partition", (part,)
    )
    return out if return_result else out.classification


def oscillation(expr, var, *, domain=True, return_result=False, partition=None):
    """Classify real oscillation on the requested domain."""
    return _function_oscillation_impl(
        expr,
        var,
        domain=domain,
        return_result=return_result,
        partition=partition,
        _domain_ready=False,
    )
