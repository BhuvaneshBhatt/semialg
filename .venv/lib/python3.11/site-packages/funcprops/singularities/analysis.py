"""Internal singularity analysis implementation."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S
from sympy.logic.boolalg import And, Or

from ..engine import evaluate_property, simplify_condition
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import (
    PropertyResult,
    Singularity,
    SingularityKind,
)
from .records import (
    _denominator_singularity_records,
    _function_branchpoint_records,
    _function_discontinuity_records,
    _function_pole_records,
    _piecewise_boundaries,
    _piecewise_discontinuous_at,
    _unique,
    branchcuts,
    definition_cuts,
    essential_singularities,
)
from .regularity import (
    _function_analytic_impl,
    _function_continuous_impl,
)


@dataclass(frozen=True)
class SingularityAnalysis:
    branch_cuts: tuple[Singularity, ...] = ()
    definition_cuts: tuple[Singularity, ...] = ()
    poles: tuple[Singularity, ...] = ()
    essential: tuple[Singularity, ...] = ()
    branch_points: tuple[Singularity, ...] = ()
    removable: tuple[Singularity, ...] = ()
    discontinuities: tuple[Singularity, ...] = ()
    domain_exterior: sp.Basic = S.false
    complete: bool = False
    unknown_heads: tuple[str, ...] = ()
    source: str = "structural-singularity-traversal"

    def all_records(self):
        return _unique(
            self.branch_cuts
            + self.definition_cuts
            + self.poles
            + self.essential
            + self.branch_points
            + self.removable
            + self.discontinuities
        )

    def condition(self):
        conds = []
        if self.domain_exterior is not S.false:
            conds.append(self.domain_exterior)
        for s in self.all_records():
            conds.append(And(s.condition, s.guard))
        return Or(*conds) if conds else S.false


def _conditioned(records, cond, *, assumptions=S.true):
    if cond is S.true and assumptions is S.true:
        return tuple(records)
    out = []
    for r in records:
        g = simplify_condition(And(r.guard, cond), assumptions=assumptions)
        if g is not S.false:
            out.append(Singularity(r.kind, r.condition, g, r.multiplicity, r.source))
    return tuple(out)


def _domain_exterior_condition(expr, var, real=True):
    """Boundary/exterior contribution of the supported natural domain."""
    try:
        from ..api import domain

        d = domain(expr, var, domain=S.Reals if real else S.Complexes)
    except SYMBOLIC_QUERY_ERRORS:
        return None
    if isinstance(d, PropertyResult):
        return None
    d = sp.sympify(d)
    if d is S.true:
        return S.false
    if d is S.false:
        return S.true
    try:
        return sp.simplify_logic(sp.Not(d))
    except SYMBOLIC_QUERY_ERRORS:
        return sp.Not(d)


def _primitive_known(expr, var):
    if expr.is_Atom:
        return True
    if expr.is_Add or expr.is_Mul or isinstance(expr, (sp.Pow, sp.Piecewise)):
        return True
    if expr.func in (
        sp.exp,
        sp.log,
        sp.sin,
        sp.cos,
        sp.tan,
        sp.cot,
        sp.sinh,
        sp.cosh,
        sp.tanh,
        sp.asin,
        sp.acos,
        sp.atan,
        sp.asinh,
        sp.acosh,
        sp.atanh,
        sp.Abs,
        sp.sign,
        sp.erf,
        sp.erfc,
        sp.LambertW,
        sp.besselj,
        sp.bessely,
        sp.gamma,
    ):
        return True
    try:
        fam = evaluate_property(expr, "LocallyAnalytic", var=var, mode="complex")
        if fam.primary is not None or fam.cases:
            return True
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return False


def _singularity_analysis_impl(
    expr, var, *, domain=S.Complexes, assumptions=S.true, ignore_unknown=False
):
    """Traverse an expression and return a typed, condition-aware singularity certificate.

    Unlike :func:`singularities`, this records whether traversal was
    complete, preserves Piecewise branch guards, includes natural-domain
    exterior conditions, and promotes poles of an inner expression to
    essential singularities when composed through a transcendental entire
    function with an essential singularity at infinity.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    complex_mode = domain in (S.Complexes, "complex", "complexes", "C", "CC")
    bc = []
    dc = []
    po = []
    es = []
    bp = []
    rm = []
    di = []
    unknown = []
    transcendental_entire = {sp.exp, sp.sin, sp.cos, sp.sinh, sp.cosh, sp.erf, sp.erfc}

    def visit(e, guard=S.true, inside_essential=False):
        if not getattr(e, "has", lambda *_: False)(var):
            return
        # Piecewise values are traversed only under their active branch guard.
        if isinstance(e, sp.Piecewise):
            covered = S.false
            for value, cond in e.args:
                active = (
                    sp.Not(covered) if cond is S.true else And(cond, sp.Not(covered))
                )
                active = simplify_condition(And(guard, active), assumptions=assumptions)
                if active is not S.false:
                    visit(value, active, inside_essential)
                covered = Or(covered, cond)
            for boundary in _piecewise_boundaries(e):
                dc.append(
                    Singularity(
                        SingularityKind.DEFINITION_CUT,
                        boundary,
                        guard,
                        source="piecewise-boundary",
                    )
                )
                verdict = _piecewise_discontinuous_at(e, var, boundary)
                if verdict is not False:
                    di.append(
                        Singularity(
                            SingularityKind.DISCONTINUITY,
                            boundary,
                            guard,
                            source="piecewise-boundary",
                        )
                    )
            return

        # Collect local information at this node, but condition it by traversal guard.
        if e.is_Function or isinstance(e, sp.Pow):
            if complex_mode:
                bc.extend(
                    _conditioned(
                        branchcuts(e, var, assumptions=assumptions, mode="records"),
                        guard,
                        assumptions=assumptions,
                    )
                )
                bp.extend(
                    _conditioned(
                        _function_branchpoint_records(e, var, assumptions=assumptions),
                        guard,
                        assumptions=assumptions,
                    )
                )
            dc.extend(
                _conditioned(
                    definition_cuts(
                        e,
                        var,
                        mode="complex" if complex_mode else "real",
                        assumptions=assumptions,
                    ),
                    guard,
                    assumptions=assumptions,
                )
            )
            denominator_poles, denominator_removable = _denominator_singularity_records(
                e, var
            )
            localp = _conditioned(
                _function_pole_records(
                    e,
                    var,
                    assumptions=assumptions,
                    denominator_poles=denominator_poles,
                ),
                guard,
                assumptions=assumptions,
            )
            if inside_essential:
                es.extend(
                    Singularity(
                        SingularityKind.ESSENTIAL,
                        p.condition,
                        p.guard,
                        p.multiplicity,
                        "outer-essential-at-infinity",
                    )
                    for p in localp
                )
            else:
                po.extend(localp)
            es.extend(
                _conditioned(
                    essential_singularities(e, var, assumptions=assumptions),
                    guard,
                    assumptions=assumptions,
                )
            )
            rm.extend(
                _conditioned(denominator_removable, guard, assumptions=assumptions)
            )
            if not complex_mode:
                di.extend(
                    _conditioned(
                        _function_discontinuity_records(
                            e, var, assumptions=assumptions
                        ),
                        guard,
                        assumptions=assumptions,
                    )
                )
            if not _primitive_known(e, var):
                unknown.append(getattr(e.func, "__name__", str(e.func)))
                if not ignore_unknown:
                    return

        outer_essential = inside_essential or (
            getattr(e, "func", None) in transcendental_entire
        )
        args = (
            [v for v, _ in e.args]
            if isinstance(e, sp.Piecewise)
            else getattr(e, "args", ())
        )
        for a in args:
            if getattr(a, "has", lambda *_: False)(var):
                visit(a, guard, outer_essential)

    visit(expr)
    dex = _domain_exterior_condition(expr, var, not complex_mode)
    complete = dex is not None and (ignore_unknown or not unknown)
    return SingularityAnalysis(
        branch_cuts=_unique(bc),
        definition_cuts=_unique(dc),
        poles=_unique(po),
        essential=_unique(es),
        branch_points=_unique(bp),
        removable=_unique(rm),
        discontinuities=_unique(di),
        domain_exterior=S.false if dex is None else dex,
        complete=complete,
        unknown_heads=tuple(sorted(set(unknown))),
    )


def singularity_analysis(
    expr, var, *, domain=S.Complexes, assumptions=S.true, ignore_unknown=False
):
    """Cached aggregate singularity analysis with assumptions in the cache key."""
    from ..proof_cache import property_proof_context

    return property_proof_context(expr, var, assumptions=assumptions).singularities(
        domain=domain, ignore_unknown=ignore_unknown
    )


def analytic(expr, var, *, domain=S.Reals, assumptions=S.true, return_result=False):
    """Classify analyticity on the requested domain."""
    return _function_analytic_impl(
        expr,
        var,
        domain=domain,
        assumptions=assumptions,
        return_result=return_result,
        _domain_ready=False,
    )


def continuous(expr, var, *, domain=S.Reals, assumptions=S.true, return_result=False):
    """Classify continuity on the requested domain."""
    return _function_continuous_impl(
        expr,
        var,
        domain=domain,
        assumptions=assumptions,
        return_result=return_result,
        _domain_ready=False,
    )
