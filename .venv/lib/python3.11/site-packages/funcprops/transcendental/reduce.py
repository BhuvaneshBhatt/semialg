from __future__ import annotations

import itertools
from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..real_analysis.common import condition_to_set, set_to_condition
from .periodic import PeriodicCondition, PeriodicSet, common_positive_period
from .roots import transcendental_roots

_REL = (sp.Eq, sp.Ne, sp.Lt, sp.Le, sp.Gt, sp.Ge)


@dataclass(frozen=True)
class UnivariateReductionResult:
    condition: object
    complete: bool
    cells: tuple[sp.Set, ...] = ()
    roots: tuple[sp.Expr, ...] = ()
    method: str = "univariate_cells"
    evidence: tuple[object, ...] = ()

    @property
    def certified(self):
        return self.complete


def _normalize_relations(condition):
    condition = sp.sympify(condition)
    replacements = {}
    for rel in condition.atoms(sp.core.relational.Relational):
        lhs = sp.simplify(rel.lhs - rel.rhs)
        replacements[rel] = rel.func(lhs, S.Zero)
    return condition.xreplace(replacements)


def _lhs_expressions(condition):
    return tuple(
        sorted(
            {
                sp.simplify(rel.lhs)
                for rel in condition.atoms(sp.core.relational.Relational)
            },
            key=sp.default_sort_key,
        )
    )


def _bool_at(condition, var, value):
    try:
        q = sp.simplify(condition.subs(var, value))
        if q in (S.true, S.false):
            return bool(q)
        q = sp.refine(q)
        if q in (S.true, S.false):
            return bool(q)
        return bool(q) if q.is_Boolean and q in (S.true, S.false) else None
    except SYMBOLIC_QUERY_ERRORS:
        return None


def _sample(a, b):
    if a is -sp.oo and b is sp.oo:
        return S.Zero
    if a is -sp.oo:
        return sp.simplify(b - 1)
    if b is sp.oo:
        return sp.simplify(a + 1)
    return sp.simplify((a + b) / 2)


def _singular_points(expr, var, interval):
    """Certified finite real discontinuity/singularity points in an interval."""
    try:
        from sympy.calculus.util import continuous_domain

        cd = continuous_domain(expr, var, interval)
        missing = sp.Complement(interval, cd)
        if missing is S.EmptySet:
            return (), True
        if isinstance(missing, sp.FiniteSet):
            return tuple(missing), True
        # Complements/unions that simplify to finitely many points.
        simp = sp.simplify(missing)
        if isinstance(simp, sp.FiniteSet):
            return tuple(simp), True
        return (), False
    except SYMBOLIC_QUERY_ERRORS:
        # Rational denominators give exact poles even when continuous_domain
        # declines a more complicated transcendental expression.
        try:
            den = sp.denom(sp.together(expr))
            if den != 1:
                sol = sp.solveset(sp.Eq(den, 0), var, domain=interval)
                if isinstance(sol, sp.FiniteSet):
                    return tuple(sol), True
        except SYMBOLIC_QUERY_ERRORS:
            pass
        return (), False


def _reduce_interval(condition, var, interval, allow_numerical):
    """Reduce a Boolean condition on one interval using certified event points.

    Roots and singularities partition the interval into cells. A result is marked
    complete only when every event search is complete and every cell truth value
    is certified.
    """
    lhs = _lhs_expressions(condition)
    roots = []
    evidence = []
    for f in lhs:
        singular, sing_complete = _singular_points(f, var, interval)
        roots.extend(singular)
        evidence.append(("singularities", f, singular, sing_complete))
        if not sing_complete:
            return UnivariateReductionResult(
                condition,
                False,
                evidence=tuple(evidence),
                method="incomplete_singularities",
            )
        rr = transcendental_roots(f, var, interval, allow_numerical=allow_numerical)
        evidence.append(rr)
        if not rr.complete:
            return UnivariateReductionResult(
                condition, False, evidence=tuple(evidence), method="incomplete_roots"
            )
        if rr.periodic is not None:
            # Bounded interval: let solveset/intersection enumerate the exact points.
            try:
                s = sp.solveset(sp.Eq(f, 0), var, domain=interval)
                if isinstance(s, sp.FiniteSet):
                    roots.extend(s)
                else:
                    return UnivariateReductionResult(
                        condition,
                        False,
                        evidence=tuple(evidence),
                        method="periodic_unexpanded",
                    )
            except SYMBOLIC_QUERY_ERRORS:
                return UnivariateReductionResult(
                    condition,
                    False,
                    evidence=tuple(evidence),
                    method="periodic_unexpanded",
                )
        else:
            roots.extend(rr.roots)
    roots = tuple(sorted(set(map(sp.simplify, roots)), key=sp.default_sort_key))
    interior = [
        r
        for r in roots
        if r != interval.start
        and r != interval.end
        and sp.Contains(r, interval) is not S.false
    ]
    boundaries = [interval.start, *interior, interval.end]
    cells = []
    true_parts = []
    for a, b in itertools.pairwise(boundaries):
        if a == b:
            continue
        sample = _sample(a, b)
        val = _bool_at(condition, var, sample)
        if val is None:
            return UnivariateReductionResult(
                condition,
                False,
                tuple(cells),
                roots,
                "undecidable_sample",
                tuple(evidence),
            )
        cell = sp.Interval.open(a, b)
        cells.append(cell)
        if val:
            true_parts.append(cell)
    for r in roots:
        if sp.Contains(r, interval) is S.false:
            continue
        val = _bool_at(condition, var, r)
        if val is None:
            continue
        cell = sp.FiniteSet(r)
        cells.append(cell)
        if val:
            true_parts.append(cell)
    # Include finite closed endpoints that are not roots of a relation.
    for endpoint, is_open in (
        (interval.start, interval.left_open),
        (interval.end, interval.right_open),
    ):
        if endpoint in (-sp.oo, sp.oo) or is_open or endpoint in roots:
            continue
        val = _bool_at(condition, var, endpoint)
        if val:
            true_parts.append(sp.FiniteSet(endpoint))
    solution = sp.Union(*true_parts) if true_parts else S.EmptySet
    return UnivariateReductionResult(
        set_to_condition(solution, var),
        True,
        tuple(cells),
        roots,
        "univariate_cells",
        tuple(evidence),
    )


def reduce_univariate_condition(
    condition, var, domain=S.Reals, *, allow_numerical=False
):
    """Certified real-univariate reduction using exact real-univariate reduction.

    Every relational boundary and discontinuity supplied by the caller should be
    represented as a relation in ``condition``.  The reducer isolates all roots,
    partitions the domain, samples each open cell and each boundary point, then
    reconstructs the exact solution.  It refuses certification if any root set
    is incomplete.
    """
    condition = _normalize_relations(condition)
    var = sp.sympify(var)
    dset = domain if isinstance(domain, sp.Set) else condition_to_set(domain, var)
    if dset is None:
        return UnivariateReductionResult(condition, False, method="unsupported_domain")

    lhs = _lhs_expressions(condition)
    if dset is S.Reals and lhs:
        p = common_positive_period(lhs, var)
        if p is not None:
            fi = sp.Interval.Ropen(0, p)
            base = _reduce_interval(condition, var, fi, allow_numerical)
            if base.complete:
                bset = condition_to_set(base.condition, var)
                if bset is not None:
                    ps = PeriodicSet(sp.Intersection(bset, fi), p, var, S.Zero)
                    return UnivariateReductionResult(
                        PeriodicCondition(ps),
                        True,
                        base.cells,
                        base.roots,
                        "periodic_cells",
                        base.evidence,
                    )

    intervals = []
    if isinstance(dset, sp.Interval):
        intervals = [dset]
    elif isinstance(dset, sp.Union) and all(
        isinstance(a, sp.Interval) for a in dset.args
    ):
        intervals = list(dset.args)
    else:
        return UnivariateReductionResult(
            condition, False, method="unsupported_domain_set"
        )
    parts = []
    allcells = []
    allroots = []
    evidence = []
    for interval in intervals:
        r = _reduce_interval(condition, var, interval, allow_numerical)
        evidence.extend(r.evidence)
        if not r.complete:
            return UnivariateReductionResult(
                condition,
                False,
                tuple(allcells),
                tuple(allroots),
                r.method,
                tuple(evidence),
            )
        parts.append(r.condition)
        allcells.extend(r.cells)
        allroots.extend(r.roots)
    out = sp.Or(*parts) if parts else S.false
    return UnivariateReductionResult(
        sp.simplify(out),
        True,
        tuple(allcells),
        tuple(allroots),
        "univariate_cells",
        tuple(evidence),
    )
