from __future__ import annotations

import sympy as sp
from sympy import S

from .. import _semialg
from ..engine import evaluate_property
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import (
    ConditionalExpression,
    Convexity,
    Monotonicity,
    PropertyResult,
    Sign,
)
from ..stratification import (
    evaluate_parameter_strata,
    normalize_cases,
    semantic_condition,
)
from ..transcendental.zero_geometry import property_interval_decomposition
from .advanced_convexity import (
    ConvexityCertificate,
    convexity_implications,
    parameterized_exp_convexity,
    parameterized_power_convexity,
    secondary_convexity_properties,
    special_function_convexity,
    strong_convexity,
)
from .common import (
    affine_real_parameter_guard,
    combine_domain_condition,
    condition_to_set,
    normalize_variables,
)
from .domain import domain as _domain
from .monotonicity import monotonicity
from .piecewise import certify_piecewise_boundaries, piecewise_decomposition
from .sign import sign
from .special_functions import parameter_guard, positive_axis_convexity_strata
from .staging import NO_STAGE_RESULT


def _from_code(value):
    try:
        code = int(str(value))
    except SYMBOLIC_QUERY_ERRORS:
        return Convexity.UNKNOWN
    return {
        1: Convexity.STRICTLY_CONVEX,
        2: Convexity.CONVEX,
        -1: Convexity.STRICTLY_CONCAVE,
        -2: Convexity.CONCAVE,
        0: Convexity.NEITHER,
    }.get(code, Convexity.UNKNOWN)


def _registry_convexity(expr, var):
    r = evaluate_property(expr, "GlobalConvexity", var=var, mode="real")
    pguard = parameter_guard(expr, var, "GlobalConvexity")
    if r.primary is not None:
        v = _from_code(r.primary)
        if pguard is S.true:
            return v
        return normalize_cases(
            (ConditionalExpression(v, pguard),),
            source="registry-convexity-stratification",
        )
    if r.cases:
        return normalize_cases(
            (
                ConditionalExpression(
                    _from_code(c.value), sp.And(semantic_condition(c.condition), pguard)
                )
                for c in r.cases
            ),
            source="registry-convexity-stratification",
        )
    return Convexity.UNKNOWN


def _normalize_semialg(value):
    return _semialg.normalize_convexity(value)


def _sum_closure(expr, var, domain, strict):
    classes = [
        convexity(a, var, domain=domain, strict_inequalities=strict)
        for a in expr.args
        if a.has(var)
    ]
    if not classes:
        return Convexity.AFFINE
    cv = (Convexity.STRICTLY_CONVEX, Convexity.CONVEX, Convexity.AFFINE)
    cc = (Convexity.STRICTLY_CONCAVE, Convexity.CONCAVE, Convexity.AFFINE)
    if all(c in cv for c in classes):
        return (
            Convexity.STRICTLY_CONVEX
            if any(c is Convexity.STRICTLY_CONVEX for c in classes)
            else Convexity.CONVEX
        )
    if all(c in cc for c in classes):
        return (
            Convexity.STRICTLY_CONCAVE
            if any(c is Convexity.STRICTLY_CONCAVE for c in classes)
            else Convexity.CONCAVE
        )
    return Convexity.UNKNOWN


def _composition(expr, var, domain, strict):
    varying = [(i, a) for i, a in enumerate(expr.args) if a.has(var)]
    if len(varying) != 1 or expr.func in (sp.Add, sp.Mul, sp.Pow, sp.Piecewise):
        return Convexity.UNKNOWN
    i, inner = varying[0]
    X = sp.Dummy("X", real=True)
    args = list(expr.args)
    args[i] = X
    outer = expr.func(*args)
    oc = _registry_convexity(outer, X)
    if outer.func is sp.exp and oc is Convexity.UNKNOWN:
        oc = Convexity.STRICTLY_CONVEX
    om = monotonicity(outer, X)
    ic = convexity(inner, var, domain=domain, strict_inequalities=strict)
    inc = om in (Monotonicity.STRICTLY_INCREASING, Monotonicity.NONDECREASING)
    dec = om in (Monotonicity.STRICTLY_DECREASING, Monotonicity.NONINCREASING)
    if (
        oc in (Convexity.STRICTLY_CONVEX, Convexity.CONVEX)
        and inc
        and ic in (Convexity.STRICTLY_CONVEX, Convexity.CONVEX, Convexity.AFFINE)
    ):
        return (
            Convexity.STRICTLY_CONVEX
            if oc is Convexity.STRICTLY_CONVEX or ic is Convexity.STRICTLY_CONVEX
            else Convexity.CONVEX
        )
    if (
        oc in (Convexity.STRICTLY_CONCAVE, Convexity.CONCAVE)
        and inc
        and ic in (Convexity.STRICTLY_CONCAVE, Convexity.CONCAVE, Convexity.AFFINE)
    ):
        return (
            Convexity.STRICTLY_CONCAVE
            if oc is Convexity.STRICTLY_CONCAVE or ic is Convexity.STRICTLY_CONCAVE
            else Convexity.CONCAVE
        )
    if (
        oc in (Convexity.STRICTLY_CONVEX, Convexity.CONVEX)
        and dec
        and ic in (Convexity.STRICTLY_CONCAVE, Convexity.CONCAVE)
    ):
        return Convexity.CONVEX
    if (
        oc in (Convexity.STRICTLY_CONCAVE, Convexity.CONCAVE)
        and dec
        and ic in (Convexity.STRICTLY_CONVEX, Convexity.CONVEX)
    ):
        return Convexity.CONCAVE
    return Convexity.UNKNOWN


def _convexity_from_second_sign(v):
    if v is Sign.POSITIVE:
        return Convexity.STRICTLY_CONVEX
    if v is Sign.NEGATIVE:
        return Convexity.STRICTLY_CONCAVE
    if v is Sign.NONNEGATIVE:
        return Convexity.CONVEX
    if v is Sign.NONPOSITIVE:
        return Convexity.CONCAVE
    if v is Sign.ZERO:
        return Convexity.AFFINE
    if v is Sign.SIGN_CHANGING:
        return Convexity.NEITHER
    return Convexity.UNKNOWN


def _flip_convexity(v):
    return {
        Convexity.STRONGLY_CONVEX: Convexity.STRONGLY_CONCAVE,
        Convexity.STRICTLY_CONVEX: Convexity.STRICTLY_CONCAVE,
        Convexity.CONVEX: Convexity.CONCAVE,
        Convexity.AFFINE: Convexity.AFFINE,
        Convexity.STRONGLY_CONCAVE: Convexity.STRONGLY_CONVEX,
        Convexity.STRICTLY_CONCAVE: Convexity.STRICTLY_CONVEX,
        Convexity.CONCAVE: Convexity.CONVEX,
    }.get(v, Convexity.UNKNOWN)


def _constant_hessian_convexity(expr, variables):
    """Exact classification for expressions with a constant symmetric Hessian.

    This intentionally handles only the zero-cost exact case. Variable
    Hessians remain on the established scalar/semialgebraic paths
    triggering expensive multivariate inequality solving here.
    """
    vv = tuple(variables)
    try:
        H = sp.hessian(expr, vv)
    except SYMBOLIC_QUERY_ERRORS:
        return Convexity.UNKNOWN
    if any(entry.has(*vv) for entry in H):
        return Convexity.UNKNOWN
    if sp.zeros(len(vv)) == H:
        return Convexity.AFFINE
    try:
        eigs = H.eigenvals()
    except SYMBOLIC_QUERY_ERRORS:
        eigs = {}
    if eigs:
        vals = list(eigs)
        if all(v.is_positive is True for v in vals):
            return Convexity.STRONGLY_CONVEX
        if all(v.is_negative is True for v in vals):
            return Convexity.STRONGLY_CONCAVE
        if all(v.is_nonnegative is True for v in vals):
            return Convexity.CONVEX
        if all(v.is_nonpositive is True for v in vals):
            return Convexity.CONCAVE
        if any(v.is_positive is True for v in vals) and any(
            v.is_negative is True for v in vals
        ):
            return Convexity.NEITHER
    # Matrix definiteness predicates are exact/tri-state; never use floating
    # eigenvalue signs to turn an undecided symbolic proof into a certificate.
    try:
        if H.is_positive_definite is True:
            return Convexity.STRONGLY_CONVEX
        if H.is_negative_definite is True:
            return Convexity.STRONGLY_CONCAVE
        if H.is_positive_semidefinite is True:
            return Convexity.CONVEX
        if H.is_negative_semidefinite is True:
            return Convexity.CONCAVE
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return Convexity.UNKNOWN


def _piecewise_convexity(expr, var, domain, strict):
    """Combine convex branch cells after certified value/slope seam checks."""
    if not isinstance(expr, sp.Piecewise):
        return Convexity.UNKNOWN
    decomp = piecewise_decomposition(expr, var, domain)
    branches = tuple(b for b in decomp.branches if isinstance(b.domain, sp.Interval))
    if not branches:
        return Convexity.UNKNOWN
    classes = [
        convexity(
            branch.expression,
            var,
            domain=branch.condition,
            strict_inequalities=strict,
        )
        for branch in branches
    ]
    if any(isinstance(value, PropertyResult) for value in classes):
        return Convexity.UNKNOWN
    convex = {
        Convexity.STRONGLY_CONVEX,
        Convexity.STRICTLY_CONVEX,
        Convexity.CONVEX,
        Convexity.AFFINE,
    }
    concave = {
        Convexity.STRONGLY_CONCAVE,
        Convexity.STRICTLY_CONCAVE,
        Convexity.CONCAVE,
        Convexity.AFFINE,
    }
    orientation = (
        Convexity.CONVEX
        if all(value in convex for value in classes)
        else Convexity.CONCAVE
        if all(value in concave for value in classes)
        else Convexity.UNKNOWN
    )
    if orientation is Convexity.UNKNOWN or not decomp.complete:
        return orientation if len(branches) == 1 else Convexity.UNKNOWN
    seams = certify_piecewise_boundaries(decomp, max_order=1)
    for seam in seams:
        if seam.continuous is False:
            return Convexity.NEITHER
        if seam.continuous is not True or not seam.derivative_limits:
            return Convexity.UNKNOWN
        left_slope, right_slope = seam.derivative_limits[0]
        if left_slope is None or right_slope is None:
            return Convexity.UNKNOWN
        jump = sp.simplify(right_slope - left_slope)
        if orientation is Convexity.CONVEX:
            if jump.is_nonnegative is not True:
                return (
                    Convexity.NEITHER if jump.is_negative is True else Convexity.UNKNOWN
                )
        elif jump.is_nonpositive is not True:
            return Convexity.NEITHER if jump.is_positive is True else Convexity.UNKNOWN
    return orientation


def _convexity_property_result(
    primary,
    domain,
    source,
    *,
    all_properties=(),
    evidence=(),
):
    """Build one normalized detailed convexity result and proof certificate."""
    proof = tuple(evidence)
    if primary is not Convexity.UNKNOWN and not any(
        isinstance(item, ConvexityCertificate) for item in proof
    ):
        proof = (
            ConvexityCertificate(primary, domain, source, proof=proof),
            *proof,
        )
    return PropertyResult(
        primary,
        all_properties=tuple(all_properties),
        source=source,
        evidence=proof,
    )


def _semialgebraic_convexity_stage(
    expr, var, assumptions, *, strict_inequalities, options
):
    """Return a decisive exact semialgebraic convexity class or ``None``."""
    if not _semialg.eligible_problem(expr, var, domain=assumptions):
        return None
    try:
        raw = _semialg.convexity(
            expr,
            var,
            domain=assumptions,
            strict_inequalities=strict_inequalities,
            **options,
        )
    except SYMBOLIC_QUERY_ERRORS:
        return None
    value = _normalize_semialg(raw)
    return value if value is not Convexity.UNKNOWN else None


def _multivariate_convexity_stage(
    expr, var, *, all_properties, return_result, intervals
):
    if not isinstance(var, (tuple, list, sp.Tuple)):
        return NO_STAGE_RESULT
    variables = normalize_variables(var)
    primary = (
        Convexity.AFFINE
        if not any(expr.has(v) for v in variables)
        else _constant_hessian_convexity(expr, variables)
    )
    detailed = all_properties or return_result or intervals
    if primary is Convexity.UNKNOWN:
        return (
            _convexity_property_result(
                primary, S.Reals, "multivariate-variable-hessian"
            )
            if detailed
            else primary
        )
    return (
        _convexity_property_result(
            primary,
            S.Reals,
            "constant-hessian",
            all_properties=convexity_implications(primary),
        )
        if detailed
        else primary
    )


def _convexity_domain_stage(assumptions, domain, var, *, detailed):
    domain_set = condition_to_set(assumptions, var)
    requested_set = (
        condition_to_set(domain, var)
        if domain not in (True, None, S.true, S.Reals)
        else None
    )
    if isinstance(requested_set, sp.Union) or isinstance(domain_set, sp.Union):
        result = Convexity.NEITHER
        return None, (
            _convexity_property_result(result, assumptions, "nonconvex-domain")
            if detailed
            else result
        )
    if domain not in (True, None, S.true, S.Reals) and not isinstance(
        requested_set, sp.Set
    ):
        result = Convexity.UNKNOWN
        return None, (
            _convexity_property_result(
                result, assumptions, "unresolved-domain-convexity"
            )
            if detailed
            else result
        )
    geometry_domain = (
        domain_set
        if isinstance(domain_set, sp.Set)
        else requested_set
        if isinstance(requested_set, sp.Set)
        else S.Reals
    )
    return geometry_domain, NO_STAGE_RESULT


def _convexity_fallback_stage(
    expr,
    var,
    *,
    domain,
    assumptions,
    strict_inequalities,
    all_properties,
    return_result,
    intervals,
    curvature_decomposition,
    options,
):
    registry_hint = Convexity.UNKNOWN
    primary = Convexity.UNKNOWN
    source = "insufficient-evidence"
    evidence = ()

    backend = _semialgebraic_convexity_stage(
        expr,
        var,
        assumptions,
        strict_inequalities=strict_inequalities,
        options=options,
    )
    if backend is not None:
        primary = backend
        source = "semialg"

    if primary is Convexity.UNKNOWN and domain is True:
        registry_hint = _registry_convexity(expr, var)
        if isinstance(registry_hint, PropertyResult):
            # Parameterized registry facts are already exact strata.
            return registry_hint
        if registry_hint not in (Convexity.UNKNOWN, Convexity.NEITHER):
            primary = registry_hint
            source = "registry"
        elif expr.func is sp.exp:
            primary = Convexity.STRICTLY_CONVEX
            source = "primitive"

    if primary is Convexity.UNKNOWN and expr.func is sp.Add:
        c = _sum_closure(expr, var, assumptions, strict_inequalities)
        if c is not Convexity.UNKNOWN:
            primary = c
            source = "sum-closure"
    if primary is Convexity.UNKNOWN and expr.func is sp.Mul:
        const, dep = expr.as_independent(var, as_Add=False)
        if const != 1 and dep != expr:
            c = _function_convexity_impl(
                dep, var, domain=assumptions, strict_inequalities=strict_inequalities
            )
            sg = sign(const, var)
            if isinstance(c, PropertyResult) or isinstance(sg, PropertyResult):
                ccases = (
                    c.cases
                    if isinstance(c, PropertyResult)
                    else (ConditionalExpression(c, S.true),)
                )
                scases = (
                    sg.cases
                    if isinstance(sg, PropertyResult)
                    else (ConditionalExpression(sg, S.true),)
                )
                cases = []
                for cv in ccases:
                    for sv in scases:
                        val = (
                            cv.value
                            if sv.value is Sign.POSITIVE
                            else _flip_convexity(cv.value)
                            if sv.value is Sign.NEGATIVE
                            else Convexity.AFFINE
                            if sv.value is Sign.ZERO
                            else Convexity.UNKNOWN
                        )
                        if val is not Convexity.UNKNOWN:
                            cases.append(
                                ConditionalExpression(
                                    val, sp.And(cv.condition, sv.condition)
                                )
                            )
                if cases:
                    return normalize_cases(
                        cases, source="scalar-convexity-stratification"
                    )
            else:
                if sg is Sign.POSITIVE:
                    primary = c
                    source = "scalar-closure"
                elif sg is Sign.NEGATIVE:
                    primary = _flip_convexity(c)
                    source = "scalar-closure"
                elif sg is Sign.ZERO:
                    primary = Convexity.AFFINE
                    source = "scalar-closure"

    if primary is Convexity.UNKNOWN:
        c = _composition(expr, var, assumptions, strict_inequalities)
        if c is not Convexity.UNKNOWN:
            primary = c
            source = "composition"

    if primary is Convexity.UNKNOWN:
        try:
            from ..proof_cache import cached_derivative

            d2 = cached_derivative(expr, var, order=2)
        except SYMBOLIC_QUERY_ERRORS:
            d2 = None
        if d2 == 0:
            primary = Convexity.AFFINE
            source = "second-derivative"
        elif d2 is not None:
            sg = sign(d2, var, domain=assumptions, strict_inequalities=True)
            if isinstance(sg, PropertyResult):
                out = normalize_cases(
                    (
                        ConditionalExpression(
                            _convexity_from_second_sign(c.value), c.condition
                        )
                        for c in sg.cases
                        if _convexity_from_second_sign(c.value) is not Convexity.UNKNOWN
                    ),
                    source="second-derivative-stratification",
                )
                if out.cases:
                    return out
            primary = _convexity_from_second_sign(sg)
            if primary is not Convexity.UNKNOWN:
                source = "second-derivative"
    if primary is Convexity.UNKNOWN and registry_hint is Convexity.NEITHER:
        primary = Convexity.NEITHER
        source = "registry"

    # Upgrade an ordinary classification when a positive curvature bound is
    # available. If an upstream backend already proved strong convexity, compute
    # the local certificate only when the caller requested proof return_result.
    ordinary = (
        Convexity.STRICTLY_CONVEX,
        Convexity.CONVEX,
        Convexity.STRICTLY_CONCAVE,
        Convexity.CONCAVE,
    )
    strong = (Convexity.STRONGLY_CONVEX, Convexity.STRONGLY_CONCAVE)
    need_certificate = (
        all_properties or return_result or intervals
    ) and primary in strong
    if primary in ordinary or need_certificate:
        strong_result = strong_convexity(expr, var, assumptions)
        if strong_result is not None:
            proved_class, certificate = strong_result
            if primary in ordinary or proved_class is primary:
                primary = proved_class
                evidence = (certificate,)
                source = "strong-convexity-bound"

    specials = special_function_convexity(expr, var, assumptions)
    if primary is Convexity.UNKNOWN:
        for q in specials:
            if q in (
                Convexity.STRONGLY_CONVEX,
                Convexity.STRICTLY_CONVEX,
                Convexity.CONVEX,
                Convexity.STRONGLY_CONCAVE,
                Convexity.STRICTLY_CONCAVE,
                Convexity.CONCAVE,
            ):
                primary = q
                source = "special-function"
                break

    if all_properties or return_result or intervals:
        extra, _ = secondary_convexity_properties(expr, var, assumptions, primary)
        for q in specials:
            if q is not primary and q not in extra:
                extra = (*extra, q)
        if (
            curvature_decomposition is not None
            and curvature_decomposition.boundary_geometry
        ):
            evidence = (*evidence, curvature_decomposition)
        return _convexity_property_result(
            primary,
            assumptions,
            source,
            all_properties=extra,
            evidence=evidence,
        )
    return primary


def _function_convexity_impl(
    expr,
    var,
    *,
    domain=True,
    strict_inequalities=True,
    all_properties=False,
    return_result=False,
    intervals=False,
    _domain_ready=False,
    **kwargs,
):
    """Classify real convexity with strong/log/quasi/pseudo refinements.

    The default is the strongest canonical ordinary convexity/concavity class.
    With ``all_properties=True`` a PropertyResult includes additional proved
    log-, quasi-, and pseudo-convexity/concavity properties.
    """
    expr = sp.sympify(expr)
    multivariate = _multivariate_convexity_stage(
        expr,
        var,
        all_properties=all_properties,
        return_result=return_result,
        intervals=intervals,
    )
    if multivariate is not NO_STAGE_RESULT:
        return multivariate
    var = sp.sympify(var)
    if expr.func is sp.exp:
        affine = affine_real_parameter_guard(expr.args[0], var)
        if affine is not None and affine[2] is S.false:
            return Convexity.UNKNOWN
    if expr.func is sp.Pow and expr.exp.is_real is False:
        return Convexity.UNKNOWN
    if not expr.has(var):
        primary = Convexity.AFFINE
        return (
            _convexity_property_result(
                primary,
                S.Reals,
                "constant",
                all_properties=(Convexity.CONVEX, Convexity.CONCAVE),
            )
            if (all_properties or return_result)
            else primary
        )
    natural = _domain(expr, var)
    assumptions = (
        sp.sympify(domain)
        if _domain_ready
        else combine_domain_condition(natural, domain, var)
    )
    if not _domain_ready:
        split = evaluate_parameter_strata(
            assumptions,
            var,
            lambda branch_domain: _function_convexity_impl(
                expr,
                var,
                domain=branch_domain,
                strict_inequalities=strict_inequalities,
                all_properties=all_properties,
                return_result=return_result,
                intervals=intervals,
                _domain_ready=True,
                **kwargs,
            ),
            source="convexity-parameter-stratification",
        )
        if split is not None:
            return split
    geom_domain, domain_result = _convexity_domain_stage(
        assumptions,
        domain,
        var,
        detailed=all_properties or return_result or intervals,
    )
    if domain_result is not NO_STAGE_RESULT:
        return domain_result
    piecewise = _piecewise_convexity(expr, var, assumptions, strict_inequalities)
    if piecewise is not Convexity.UNKNOWN:
        if all_properties or return_result or intervals:
            return _convexity_property_result(
                piecewise,
                assumptions,
                "piecewise-analysis",
                all_properties=convexity_implications(piecewise),
            )
        return piecewise

    # Cheap theorem strata precede zero-geometry/semialgebraic analysis.  This is
    # important for special functions whose derivatives are transcendental but
    # whose curvature is known exactly.
    strat = parameterized_power_convexity(
        expr, var, assumptions
    ) or parameterized_exp_convexity(expr, var)
    if strat is not None:
        return strat
    try:
        positive_axis = (
            sp.simplify_logic(sp.Implies(semantic_condition(assumptions), var > 0))
            is S.true
        )
    except SYMBOLIC_QUERY_ERRORS:
        positive_axis = False
    if positive_axis:
        sfacts = positive_axis_convexity_strata(expr, var)
        if sfacts:
            out = normalize_cases(
                (ConditionalExpression(v, g) for v, g, *_ in sfacts),
                source="special-function-convexity-stratification",
            )
            if out.primary is not None:
                if all_properties or return_result or intervals:
                    extras = tuple(
                        q
                        for q in special_function_convexity(expr, var, assumptions)
                        if q is not out.primary
                    )
                    return _convexity_property_result(
                        out.primary,
                        assumptions,
                        out.source or "special-function",
                        all_properties=extras,
                    )
                return out.primary
            return out
    curvature_decomposition = (
        property_interval_decomposition(
            expr,
            var,
            derivative_order=2,
            domain=geom_domain,
            positive_value=Convexity.STRICTLY_CONVEX,
            negative_value=Convexity.STRICTLY_CONCAVE,
        )
        if (intervals or return_result)
        else None
    )

    # Certified second-derivative zero geometry can establish alternating
    # convex/concave cells even when symbolic sign simplification cannot.
    raw_curv = property_interval_decomposition(
        expr,
        var,
        derivative_order=2,
        domain=geom_domain,
        positive_value=Convexity.STRICTLY_CONVEX,
        negative_value=Convexity.STRICTLY_CONCAVE,
    )
    if (
        raw_curv.complete
        and raw_curv.alternating
        and raw_curv.boundary_geometry
        and geom_domain is S.Reals
    ):
        if intervals or return_result or all_properties:
            return _convexity_property_result(
                Convexity.NEITHER,
                assumptions,
                "second-derivative-zero-geometry",
                evidence=(raw_curv,),
            )
        return Convexity.NEITHER

    return _convexity_fallback_stage(
        expr,
        var,
        domain=domain,
        assumptions=assumptions,
        strict_inequalities=strict_inequalities,
        all_properties=all_properties,
        return_result=return_result,
        intervals=intervals,
        curvature_decomposition=curvature_decomposition,
        options=kwargs,
    )


def convexity(
    expr,
    var,
    *,
    domain=True,
    strict_inequalities=True,
    all_properties=False,
    return_result=False,
    intervals=False,
    **kwargs,
):
    """Classify real convexity with optional detailed property metadata."""
    return _function_convexity_impl(
        expr,
        var,
        domain=domain,
        strict_inequalities=strict_inequalities,
        all_properties=all_properties,
        return_result=return_result,
        intervals=intervals,
        _domain_ready=False,
        **kwargs,
    )
