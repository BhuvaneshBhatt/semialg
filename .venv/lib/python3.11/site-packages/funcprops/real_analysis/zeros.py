from __future__ import annotations

import sympy as sp
from sympy import S

from ..exactness import certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import PropertyResult
from ..stratification import evaluate_parameter_strata
from ..transcendental.roots import (
    RootCertificate,
    RootIsolationResult,
    transcendental_roots,
)
from .common import condition_to_set
from .domain import domain as _domain


def _merge_root_results(results, domain_set):
    """Merge per-cell zero results without discarding certification metadata."""
    roots = []
    zero_intervals = []
    isolated = []
    certificates = []
    evidence = []
    families = []
    geometry = []
    complete = True
    methods = []
    periodic = None
    for result in results:
        roots.extend(result.roots)
        zero_intervals.extend(result.zero_intervals)
        isolated.extend(result.isolated_intervals)
        certificates.extend(result.certificates)
        evidence.extend(result.evidence)
        families.extend(result.parametric_families)
        geometry.extend(result.geometry)
        complete = complete and result.complete
        methods.append(result.method)
        if result.periodic is not None:
            # A periodic representation cannot in general be combined losslessly
            # with unrelated domain cells. Preserve it only for the single-cell case.
            if periodic is None and len(results) == 1:
                periodic = result.periodic
            else:
                complete = False
    roots = tuple(sorted(set(roots), key=sp.default_sort_key))
    zero_intervals = tuple(sorted(set(zero_intervals), key=sp.default_sort_key))
    method = "domain_cell_isolation:" + "+".join(dict.fromkeys(methods))
    return RootIsolationResult(
        roots,
        zero_intervals,
        tuple(isolated),
        periodic,
        complete,
        method,
        (*tuple(evidence), ("domain_cells", domain_set)),
        tuple(certificates),
        tuple(families),
        tuple(geometry),
    )


def _domain_cells(domain_set):
    """Return explicit real interval/point cells, or ``None`` if unresolved."""
    domain_set = sp.simplify(domain_set)
    if domain_set is S.EmptySet:
        return ()
    if isinstance(domain_set, (sp.Interval, sp.FiniteSet)):
        return (domain_set,)
    if isinstance(domain_set, sp.Union):
        cells = []
        for arg in domain_set.args:
            sub = _domain_cells(arg)
            if sub is None:
                return None
            cells.extend(sub)
        return tuple(cells)
    return None


def _finite_real_solutions(equation, var):
    try:
        sol = sp.solveset(equation, var, domain=S.Reals)
    except SYMBOLIC_QUERY_ERRORS:
        return ()
    if not isinstance(sol, sp.FiniteSet):
        return ()
    return tuple(q for q in sol if q.is_real is not False)


def _structural_exceptional_points(expr, var, domain_set):
    """Exact points at which a real-expression regime or branch can change."""
    equations = []
    try:
        den = sp.denom(sp.together(expr))
        if den != 1:
            equations.append(sp.Eq(den, 0))
    except SYMBOLIC_QUERY_ERRORS:
        pass
    for node in sp.preorder_traversal(expr):
        if node.func is sp.log and node.args:
            equations.append(sp.Eq(node.args[0], 0))
        elif isinstance(node, sp.Pow):
            base, exponent = node.args
            if exponent.is_Rational and not exponent.is_Integer:
                equations.append(sp.Eq(base, 0))
        elif isinstance(node, sp.Piecewise):
            for _value, cond in node.args:
                if cond is True or cond is S.true:
                    continue
                from sympy.core.relational import Relational

                for rel in cond.atoms(Relational):
                    equations.append(sp.Eq(rel.lhs - rel.rhs, 0))
    points = []
    for eq in equations:
        for point in _finite_real_solutions(eq, var):
            try:
                contained = domain_set.contains(point)
            except SYMBOLIC_QUERY_ERRORS:
                continue
            if contained is S.true and point not in points:
                points.append(point)
    return tuple(sorted(points, key=sp.default_sort_key))


def _split_interval_at_points(interval, points):
    if not isinstance(interval, sp.Interval):
        return (interval,)
    interior = []
    for point in points:
        try:
            if interval.contains(point) is S.true and point not in (
                interval.start,
                interval.end,
            ):
                interior.append(point)
        except SYMBOLIC_QUERY_ERRORS:
            pass
    interior = tuple(sorted(set(interior), key=sp.default_sort_key))
    if not interior:
        return (interval,)
    cells = []
    left = interval.start
    left_open = interval.left_open
    for point in interior:
        cells.append(sp.Interval(left, point, left_open=left_open, right_open=True))
        cells.append(sp.FiniteSet(point))
        left = point
        left_open = True
    cells.append(
        sp.Interval(
            left, interval.end, left_open=left_open, right_open=interval.right_open
        )
    )
    return tuple(c for c in cells if c is not S.EmptySet)


def _zeros_on_set(expr, var, domain_set, *, return_result, allow_numerical):
    """Evaluate zero geometry on one parameter-free real set.

    Disconnected real domains are analyzed cell by cell so no certificate spans a
    pole, branch boundary, logarithmic exclusion, or other point where the input is
    undefined.
    """
    cells = _domain_cells(domain_set)
    if cells is not None:
        exceptional = _structural_exceptional_points(expr, var, domain_set)
        split_cells = []
        for cell in cells:
            split_cells.extend(_split_interval_at_points(cell, exceptional))
        cells = tuple(split_cells)
    if cells is not None and not (
        len(cells) == 1 and isinstance(cells[0], sp.Interval) and cells[0] == domain_set
    ):
        results = []
        for cell in cells:
            if isinstance(cell, sp.FiniteSet):
                roots = []
                certs = []
                cell_complete = True
                for point in cell:
                    try:
                        value = sp.simplify(expr.subs(var, point))
                    except SYMBOLIC_QUERY_ERRORS:
                        cell_complete = False
                        continue
                    zero = certify_zero(value)
                    if zero is True:
                        roots.append(point)
                        certs.append(
                            RootCertificate(
                                sp.FiniteSet(point),
                                1,
                                "exact_domain_point_root",
                                (("value", value),),
                                True,
                            )
                        )
                    elif zero is None:
                        cell_complete = False
                results.append(
                    RootIsolationResult(
                        tuple(roots),
                        complete=cell_complete,
                        method="domain_point_evaluation",
                        certificates=tuple(certs),
                    )
                )
            else:
                results.append(
                    transcendental_roots(
                        expr, var, cell, allow_numerical=allow_numerical
                    )
                )
        rr = _merge_root_results(results, domain_set)
    else:
        rr = transcendental_roots(
            expr, var, domain_set, allow_numerical=allow_numerical
        )
    if return_result:
        return rr
    if rr.complete:
        if rr.periodic is not None:
            return rr.periodic
        if (
            rr.method == "dlmf_indexed_zeros"
            and rr.evidence
            and isinstance(rr.evidence[0], sp.Set)
        ):
            return rr.evidence[0]
        if rr.isolated_intervals or rr.parametric_families:
            return rr
        parts = []
        if rr.roots:
            parts.append(sp.FiniteSet(*rr.roots))
        parts.extend(rr.zero_intervals)
        if not parts:
            return S.EmptySet
        return sp.Union(*parts) if len(parts) > 1 else parts[0]
    return rr


def _piecewise_zero_set(expr, var, domain_set, *, return_result, allow_numerical):
    from .piecewise import piecewise_decomposition

    decomp = piecewise_decomposition(expr, var, domain_set)
    if not decomp.complete:
        return None
    branch_results = []
    plain_sets = []
    for branch in decomp.branches:
        if branch.domain is None:
            return None
        value = zeros(
            branch.expression,
            var,
            domain=branch.domain,
            return_result=return_result,
            allow_numerical=allow_numerical,
        )
        if return_result:
            if not isinstance(value, RootIsolationResult):
                return None
            branch_results.append(value)
        else:
            if not isinstance(value, sp.Set):
                return None
            plain_sets.append(value)
    if return_result:
        return _merge_root_results(tuple(branch_results), domain_set)
    return sp.Union(*plain_sets) if plain_sets else S.EmptySet


def _function_zeros_impl(
    expr,
    var,
    *,
    domain=True,
    return_result=False,
    allow_numerical=False,
    _domain_ready=False,
):
    """Return the certified real zero set of ``expr``.

    Exact finite and periodic families are returned directly. Parameter-dependent
    domain conditions are evaluated stratum by stratum. If a root is only
    certified by an isolating interval, ``RootIsolationResult`` is returned so
    the certificate is not discarded.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    natural = _domain(expr, var)
    if isinstance(natural, PropertyResult):
        # Unknown natural-domain information must not erase an explicit user
        # restriction. Fall back to the requested set/condition only.
        if domain is True:
            domain_set = S.Reals
        elif isinstance(domain, sp.Set):
            domain_set = domain
        else:
            domain_set = condition_to_set(sp.sympify(domain), var)
            if domain_set is None:
                return RootIsolationResult(
                    (), complete=False, method="unsupported-domain"
                )
        return _zeros_on_set(
            expr,
            var,
            domain_set,
            return_result=return_result,
            allow_numerical=allow_numerical,
        )

    if domain is True:
        condition = natural
    elif isinstance(domain, sp.Set):
        natural_set = condition_to_set(natural, var)
        domain_set = sp.Intersection(
            natural_set if natural_set is not None else S.Reals, domain
        )
        return _zeros_on_set(
            expr,
            var,
            domain_set,
            return_result=return_result,
            allow_numerical=allow_numerical,
        )
    else:
        condition = sp.And(natural, sp.sympify(domain))

    if not _domain_ready:
        stratified = evaluate_parameter_strata(
            condition,
            var,
            lambda branch: _function_zeros_impl(
                expr,
                var,
                domain=branch,
                return_result=return_result,
                allow_numerical=allow_numerical,
                _domain_ready=True,
            ),
            source="zero-parameter-stratification",
        )
        if stratified is not None:
            return stratified

    domain_set = condition_to_set(condition, var)
    if domain_set is None:
        return RootIsolationResult((), complete=False, method="unsupported-domain")
    if isinstance(expr, sp.Piecewise):
        merged = _piecewise_zero_set(
            expr,
            var,
            domain_set,
            return_result=return_result,
            allow_numerical=allow_numerical,
        )
        if merged is not None:
            return merged
    return _zeros_on_set(
        expr,
        var,
        domain_set,
        return_result=return_result,
        allow_numerical=allow_numerical,
    )


def zeros(expr, var, *, domain=True, return_result=False, allow_numerical=False):
    """Return the certified real zeros of ``expr`` on the requested domain."""
    return _function_zeros_impl(
        expr,
        var,
        domain=domain,
        return_result=return_result,
        allow_numerical=allow_numerical,
        _domain_ready=False,
    )
