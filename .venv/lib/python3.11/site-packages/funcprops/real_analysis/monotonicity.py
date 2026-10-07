from __future__ import annotations

import contextlib

import sympy as sp
from sympy import S

from .. import _semialg
from ..engine import evaluate_property
from ..exactness import certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import ConditionalExpression, Monotonicity, PropertyResult, Sign
from ..stratification import (
    evaluate_parameter_strata,
    normalize_cases,
    semantic_condition,
)
from .common import (
    affine_real_parameter_guard,
    combine_domain_condition,
    condition_to_set,
    is_full_natural_domain,
)
from .domain import domain as _domain
from .monotonic_certificates import (
    monotonic_decomposition,
    monotonic_partition,
    quick_nonmonotonicity_witness,
)
from .piecewise import (
    certify_piecewise_boundaries,
    combine_piecewise_monotonicity,
    piecewise_decomposition,
)
from .sign import sign
from .special_functions import parameter_guard, positive_axis_monotonicity_strata
from .staging import NO_STAGE_RESULT


def _from_code(value):
    try:
        code = int(str(value))
    except SYMBOLIC_QUERY_ERRORS:
        return Monotonicity.UNKNOWN
    return {
        1: Monotonicity.STRICTLY_INCREASING,
        2: Monotonicity.NONDECREASING,
        -1: Monotonicity.STRICTLY_DECREASING,
        -2: Monotonicity.NONINCREASING,
        0: Monotonicity.NOT_MONOTONE,
    }.get(code, Monotonicity.UNKNOWN)


def _direct_registry(expr, var):
    r = evaluate_property(expr, "GlobalMonotonicity", var=var, mode="real")
    pguard = parameter_guard(expr, var, "GlobalMonotonicity")
    if r.primary is not None:
        v = _from_code(r.primary)
        if pguard is S.true:
            return v
        return normalize_cases(
            (ConditionalExpression(v, pguard),),
            source="registry-monotonicity-stratification",
        )
    if r.cases:
        out = normalize_cases(
            (
                ConditionalExpression(
                    _from_code(c.value), sp.And(semantic_condition(c.condition), pguard)
                )
                for c in r.cases
            ),
            source="registry-monotonicity-stratification",
        )
        if out.cases:
            return out
    return Monotonicity.UNKNOWN


def _compose_monotonicity(expr, var, domain):
    # g(h(x)) with one varying argument: combine outer and inner order.
    varying = [
        (i, a)
        for i, a in enumerate(getattr(expr, "args", ()))
        if getattr(a, "has", lambda *_: False)(var)
    ]
    if len(varying) != 1 or expr.func in (sp.Add, sp.Mul, sp.Pow, sp.Piecewise):
        return Monotonicity.UNKNOWN
    i, inner = varying[0]
    # A primitive f(x) is not a composition f(g(x)); treating it as f∘id can
    # re-import an overbroad global catalogue fact on a restricted domain.
    if inner == var:
        return Monotonicity.UNKNOWN
    X = sp.Dummy("X", real=True)
    args = list(expr.args)
    args[i] = X
    outer = expr.func(*args)
    om = _direct_registry(outer, X)
    if om is Monotonicity.UNKNOWN and outer.func is sp.exp:
        om = Monotonicity.STRICTLY_INCREASING
    im = monotonicity(inner, var, domain=domain)
    if isinstance(om, PropertyResult) or isinstance(im, PropertyResult):
        return Monotonicity.UNKNOWN
    inc = (Monotonicity.STRICTLY_INCREASING, Monotonicity.NONDECREASING)
    dec = (Monotonicity.STRICTLY_DECREASING, Monotonicity.NONINCREASING)
    if om in inc and im in inc:
        return (
            Monotonicity.STRICTLY_INCREASING
            if om is inc[0] and im is inc[0]
            else Monotonicity.NONDECREASING
        )
    if om in inc and im in dec:
        return (
            Monotonicity.STRICTLY_DECREASING
            if om is inc[0] and im is dec[0]
            else Monotonicity.NONINCREASING
        )
    if om in dec and im in inc:
        return (
            Monotonicity.STRICTLY_DECREASING
            if om is dec[0] and im is inc[0]
            else Monotonicity.NONINCREASING
        )
    if om in dec and im in dec:
        return (
            Monotonicity.STRICTLY_INCREASING
            if om is dec[0] and im is dec[0]
            else Monotonicity.NONDECREASING
        )
    return Monotonicity.UNKNOWN


def _closure_monotonicity(expr, var, domain):
    """Cheap sign-aware closure rules for sums, products and Min/Max."""
    inc = (Monotonicity.STRICTLY_INCREASING, Monotonicity.NONDECREASING)
    dec = (Monotonicity.STRICTLY_DECREASING, Monotonicity.NONINCREASING)
    if expr.is_Add:
        parts = [a for a in expr.args if a.has(var)]
        if not parts:
            return Monotonicity.CONSTANT
        vals = [monotonicity(a, var, domain=domain) for a in parts]
        if any(isinstance(v, PropertyResult) for v in vals):
            return Monotonicity.UNKNOWN
        if all(v in inc for v in vals):
            return (
                Monotonicity.STRICTLY_INCREASING
                if any(v is inc[0] for v in vals)
                else Monotonicity.NONDECREASING
            )
        if all(v in dec for v in vals):
            return (
                Monotonicity.STRICTLY_DECREASING
                if any(v is dec[0] for v in vals)
                else Monotonicity.NONINCREASING
            )
    if expr.func in (sp.Min, sp.Max):
        vals = [monotonicity(a, var, domain=domain) for a in expr.args if a.has(var)]
        if vals and not any(isinstance(v, PropertyResult) for v in vals):
            if all(v in inc for v in vals):
                return Monotonicity.NONDECREASING
            if all(v in dec for v in vals):
                return Monotonicity.NONINCREASING
    if expr.is_Mul:
        const, dep = expr.as_independent(var, as_Add=False)
        parts = list(sp.Mul.make_args(dep))
        if certify_zero(const) is True:
            return Monotonicity.CONSTANT
        # Product of positive factors moving in the same direction.
        if parts:
            vals = []
            signs = []
            for a in parts:
                vals.append(monotonicity(a, var, domain=domain))
                signs.append(sign(a, var, domain=domain))
            if not any(isinstance(v, PropertyResult) for v in vals):
                nonneg = all(q in (Sign.POSITIVE, Sign.NONNEGATIVE) for q in signs)
                strictly_pos = all(q is Sign.POSITIVE for q in signs)
                if nonneg and all(v in inc for v in vals):
                    # A merely nonnegative factor can vanish on an interval;
                    # then the product may have a plateau even if another
                    # factor is strictly increasing.  Strictness is safe only
                    # when every factor is strictly positive.
                    out = (
                        Monotonicity.STRICTLY_INCREASING
                        if strictly_pos and any(v is inc[0] for v in vals)
                        else Monotonicity.NONDECREASING
                    )
                    if const.is_negative is True:
                        return (
                            Monotonicity.STRICTLY_DECREASING
                            if out is Monotonicity.STRICTLY_INCREASING
                            else Monotonicity.NONINCREASING
                        )
                    if const.is_positive is True:
                        return out
                if nonneg and all(v in dec for v in vals):
                    out = (
                        Monotonicity.STRICTLY_DECREASING
                        if strictly_pos and any(v is dec[0] for v in vals)
                        else Monotonicity.NONINCREASING
                    )
                    if const.is_negative is True:
                        return (
                            Monotonicity.STRICTLY_INCREASING
                            if out is Monotonicity.STRICTLY_DECREASING
                            else Monotonicity.NONDECREASING
                        )
                    if const.is_positive is True:
                        return out
    return Monotonicity.UNKNOWN


def _piecewise_monotonicity(expr, var, domain):
    """Classify effective branches and certify ordering through every seam."""
    if not isinstance(expr, sp.Piecewise):
        return Monotonicity.UNKNOWN
    decomp = piecewise_decomposition(expr, var, domain)
    branches = decomp.branches
    if not branches:
        return Monotonicity.UNKNOWN
    values = [
        monotonicity(b.expression, var, domain=b.condition)
        for b in branches
        if isinstance(b.domain, sp.Interval)
    ]
    if not values:
        return Monotonicity.CONSTANT if len(branches) == 1 else Monotonicity.UNKNOWN
    direction = combine_piecewise_monotonicity(values)
    if direction in (Monotonicity.UNKNOWN, Monotonicity.NOT_MONOTONE):
        return direction
    if not decomp.complete:
        return direction if len(branches) == 1 else Monotonicity.UNKNOWN
    seams = certify_piecewise_boundaries(decomp, max_order=0)
    for seam in seams:
        if seam.left_limit is None or seam.right_limit is None or seam.value is None:
            return Monotonicity.UNKNOWN
        lv = sp.simplify(seam.value - seam.left_limit)
        rv = sp.simplify(seam.right_limit - seam.value)
        if direction in (Monotonicity.NONDECREASING, Monotonicity.STRICTLY_INCREASING):
            if lv.is_nonnegative is not True or rv.is_nonnegative is not True:
                if lv.is_negative is True or rv.is_negative is True:
                    return Monotonicity.NOT_MONOTONE
                return Monotonicity.UNKNOWN
        elif direction in (
            Monotonicity.NONINCREASING,
            Monotonicity.STRICTLY_DECREASING,
        ):
            if lv.is_nonpositive is not True or rv.is_nonpositive is not True:
                if lv.is_positive is True or rv.is_positive is True:
                    return Monotonicity.NOT_MONOTONE
                return Monotonicity.UNKNOWN
    return direction


def _parameter_strata(expr, var, domain):
    if expr.func is sp.Pow and expr.base == var:
        a = expr.exp
        if a.free_symbols and domain == (var > 0):
            real_a = S.true if a.is_real is True else sp.Q.real(a)
            return normalize_cases(
                (
                    ConditionalExpression(
                        Monotonicity.STRICTLY_INCREASING, sp.And(real_a, a > 0)
                    ),
                    ConditionalExpression(Monotonicity.CONSTANT, sp.Eq(a, 0)),
                    ConditionalExpression(
                        Monotonicity.STRICTLY_DECREASING, sp.And(real_a, a < 0)
                    ),
                ),
                source="power-parameter-stratification",
            )
    if expr.func is sp.exp:
        affine = affine_real_parameter_guard(expr.args[0], var)
        if affine is not None:
            a, _b, real_guard = affine
            if real_guard is S.false:
                return None
            if a.free_symbols or real_guard is not S.true:
                return normalize_cases(
                    (
                        ConditionalExpression(
                            Monotonicity.STRICTLY_INCREASING, sp.And(real_guard, a > 0)
                        ),
                        ConditionalExpression(
                            Monotonicity.CONSTANT, sp.And(real_guard, sp.Eq(a, 0))
                        ),
                        ConditionalExpression(
                            Monotonicity.STRICTLY_DECREASING, sp.And(real_guard, a < 0)
                        ),
                    ),
                    source="exp-affine-parameter-stratification",
                )
    return None


def _semialgebraic_monotonicity_stage(
    expr, var, assumptions, *, return_result, intervals
):
    """Run the exact semialgebraic monotonicity stage, if eligible.

    Interval requests deliberately continue to funcprops' shared partition
    engine because callers rely on its boundary/zero-geometry certificates.
    """
    if intervals or not _semialg.eligible_problem(expr, var, domain=assumptions):
        return None
    try:
        raw = _semialg.monotonicity(
            expr,
            var,
            domain=assumptions,
            return_result=bool(intervals or return_result),
        )
    except (_semialg.SemialgUnavailable, *SYMBOLIC_QUERY_ERRORS):
        return None
    value = _semialg.normalize_monotonicity(raw)
    if value is Monotonicity.UNKNOWN:
        return None
    classification = getattr(raw, "classification", raw)
    evidence = (raw,) if raw is not classification else ()
    return value, evidence


def _monotonicity_output(value, *, intervals, return_result, source=None, evidence=()):
    if intervals or return_result:
        return PropertyResult(
            value, source=source or "monotonicity", evidence=tuple(evidence)
        )
    return value


def _monotonicity_theorem_stage(
    expr,
    var,
    *,
    assumptions,
    domain,
    natural,
    domain_set,
    geometry_domain,
    intervals,
    return_result,
    partition,
):
    strata = _parameter_strata(expr, var, assumptions)
    if strata is not None:
        return strata

    piecewise = _piecewise_monotonicity(expr, var, assumptions)
    if piecewise is not Monotonicity.UNKNOWN:
        return _monotonicity_output(
            piecewise,
            intervals=intervals,
            return_result=return_result,
            source="piecewise-analysis",
        )

    try:
        positive_axis = (
            isinstance(domain_set, sp.Interval)
            and domain_set.start == 0
            and domain_set.left_open
            and domain_set.end is S.Infinity
        ) or sp.simplify_logic(
            sp.Implies(semantic_condition(assumptions), var > 0)
        ) is S.true
    except SYMBOLIC_QUERY_ERRORS:
        positive_axis = False
    if positive_axis:
        facts = positive_axis_monotonicity_strata(expr, var, assumptions)
        if facts:
            result = normalize_cases(
                (ConditionalExpression(value, guard) for value, guard, *_ in facts),
                source="special-function-monotonicity-stratification",
            )
            if result.primary is None:
                return result
            return _monotonicity_output(
                result.primary,
                intervals=intervals,
                return_result=return_result,
                source="special-function-monotonicity-stratification",
            )

    if domain is True or is_full_natural_domain(assumptions, natural):
        result = _direct_registry(expr, var)
        if result is not Monotonicity.UNKNOWN:
            if intervals or return_result:
                part = (
                    partition
                    if partition is not None
                    and getattr(partition, "domain", None) == geometry_domain
                    else monotonic_partition(
                        expr, var, domain=geometry_domain, return_result=True
                    )
                )
                evidence = (
                    (part.property_decomposition,)
                    if part.property_decomposition is not None
                    else ()
                ) + ((part,) if part.complete else ())
                return _monotonicity_output(
                    result,
                    intervals=intervals,
                    return_result=return_result,
                    source="registry-monotonicity",
                    evidence=evidence,
                )
            return result

    result = _closure_monotonicity(expr, var, assumptions)
    if result is not Monotonicity.UNKNOWN:
        return _monotonicity_output(
            result,
            intervals=intervals,
            return_result=return_result,
            source="closure-properties",
        )

    backend = _semialgebraic_monotonicity_stage(
        expr, var, assumptions, return_result=return_result, intervals=intervals
    )
    if backend is not None:
        value, evidence = backend
        return _monotonicity_output(
            value,
            intervals=intervals,
            return_result=return_result,
            source="semialg",
            evidence=evidence,
        )
    return NO_STAGE_RESULT


def _monotonicity_fallback_stage(
    expr, var, *, assumptions, geometry_domain, intervals, return_result, partition
):
    def finish(value, source=None, evidence=()):
        return _monotonicity_output(
            value,
            intervals=intervals,
            return_result=return_result,
            source=source,
            evidence=evidence,
        )

    dec = monotonic_decomposition(
        expr, var, strict=True, continuous=True, return_result=True
    )
    if dec is not None and dec.inner != expr and dec.inner.has(var):
        try:
            valid = sp.simplify_logic(sp.Implies(assumptions, dec.condition))
            if valid is S.true:
                im = monotonicity(dec.inner, var, domain=assumptions)
            else:
                im = PropertyResult(None)
            if not isinstance(im, PropertyResult):
                inc = (Monotonicity.STRICTLY_INCREASING, Monotonicity.NONDECREASING)
                decm = (Monotonicity.STRICTLY_DECREASING, Monotonicity.NONINCREASING)
                if im in inc + decm:
                    increasing = im in inc
                    strict = im in (inc[0], decm[0]) and dec.strict
                    if dec.direction < 0:
                        increasing = not increasing
                    out = (
                        Monotonicity.STRICTLY_INCREASING
                        if increasing and strict
                        else Monotonicity.NONDECREASING
                        if increasing
                        else Monotonicity.STRICTLY_DECREASING
                        if strict
                        else Monotonicity.NONINCREASING
                    )
                    return finish(out, "monotonic-decomposition", (dec,))
        except SYMBOLIC_QUERY_ERRORS:
            pass

    m = _compose_monotonicity(expr, var, assumptions)
    if m is not Monotonicity.UNKNOWN:
        return finish(m, "composition")

    def from_sign(v):
        if v is Sign.POSITIVE:
            return Monotonicity.STRICTLY_INCREASING
        if v is Sign.NEGATIVE:
            return Monotonicity.STRICTLY_DECREASING
        if v is Sign.NONNEGATIVE:
            return Monotonicity.NONDECREASING
        if v is Sign.NONPOSITIVE:
            return Monotonicity.NONINCREASING
        if v is Sign.ZERO:
            return Monotonicity.CONSTANT
        if v is Sign.SIGN_CHANGING:
            return Monotonicity.NOT_MONOTONE
        return Monotonicity.UNKNOWN

    def from_partition(part):
        ev = (
            (part.property_decomposition,)
            if part.property_decomposition is not None
            else ()
        ) + (part,)
        if part.indexed is not None and part.indexed.boundary_geometry:
            # Completeness of the entire requested domain is not needed to
            # disprove monotonicity: an infinite alternating critical family
            # contained inside the domain is already a certificate.
            try:
                from ..transcendental.zero_geometry import (
                    geometry_support_inside_domain,
                )

                if any(
                    geometry_support_inside_domain(g, part.domain)
                    for g in part.indexed.boundary_geometry
                ):
                    return finish(
                        Monotonicity.NOT_MONOTONE, "indexed-monotonic-partition", ev
                    )
            except SYMBOLIC_QUERY_ERRORS:
                pass
        if not part.complete:
            return None
        if part.cells:
            dirs = {c.direction for c in part.cells}
            if dirs == {Monotonicity.STRICTLY_INCREASING}:
                return finish(
                    Monotonicity.STRICTLY_INCREASING, "monotonic-partition", ev
                )
            if dirs == {Monotonicity.STRICTLY_DECREASING}:
                return finish(
                    Monotonicity.STRICTLY_DECREASING, "monotonic-partition", ev
                )
            if len(dirs) > 1:
                return finish(Monotonicity.NOT_MONOTONE, "monotonic-partition", ev)
        return None

    # Indexed zero geometry is cheap to recognize and gives stronger evidence
    # than a generic derivative-sign attempt. Details/interval requests also
    # explicitly ask for the partition certificate.
    need_part = bool(intervals or return_result)
    if not need_part:
        try:
            from ..transcendental.zero_geometry import zero_geometry

            need_part = bool(
                zero_geometry(expr, var, derivative_order=1, domain=geometry_domain)
            )
        except SYMBOLIC_QUERY_ERRORS:
            need_part = False
    part = (
        partition
        if partition is not None
        and getattr(partition, "domain", None) == geometry_domain
        else None
    )
    if need_part and part is None:
        part = monotonic_partition(
            expr, var, domain=geometry_domain, return_result=True
        )
        ans = from_partition(part)
        if ans is not None:
            return ans

    # For generic expressions, derivative sign is often cheaper than solving
    # all critical points. Strong signs finish immediately. Weak signs are
    # retained as a fallback but a partition may upgrade them to strictness.
    try:
        from ..proof_cache import cached_derivative

        derivative = part.derivative if part is not None else None
        registry_derivative_evidence = None
        if derivative is None:
            from ..registry_derivatives import registry_derivative

            cert = registry_derivative(expr, var)
            try:
                applicable = (
                    cert.complete
                    and sp.simplify_logic(
                        sp.Implies(semantic_condition(assumptions), cert.conditions)
                    )
                    is S.true
                )
            except SYMBOLIC_QUERY_ERRORS:
                applicable = False
            if applicable:
                derivative = cert.derivative
                registry_derivative_evidence = cert
            else:
                derivative = cached_derivative(expr, var, order=1)
    except SYMBOLIC_QUERY_ERRORS:
        return finish(
            Monotonicity.UNKNOWN,
            "unknown",
            (part,) if part is not None and (return_result or intervals) else (),
        )
    if derivative == 0:
        return finish(
            Monotonicity.CONSTANT,
            "zero-derivative",
            (part,) if part is not None else (),
        )
    s = sign(derivative, var, domain=assumptions, strict_inequalities=True)
    if isinstance(s, PropertyResult):
        out = normalize_cases(
            (ConditionalExpression(from_sign(c.value), c.condition) for c in s.cases),
            source="derivative-sign-stratification",
        )
        return out if out.cases else Monotonicity.UNKNOWN
    direct = from_sign(s)
    if direct in (
        Monotonicity.STRICTLY_INCREASING,
        Monotonicity.STRICTLY_DECREASING,
        Monotonicity.CONSTANT,
        Monotonicity.NOT_MONOTONE,
    ):
        return finish(
            direct,
            "derivative-sign",
            tuple(
                e
                for e in (part, registry_derivative_evidence)
                if e is not None and (return_result or intervals)
            ),
        )

    if part is None:
        part = monotonic_partition(
            expr, var, domain=geometry_domain, return_result=True
        )
        ans = from_partition(part)
        if ans is not None:
            return ans

    witness = quick_nonmonotonicity_witness(expr, var, geometry_domain)
    if witness is not None:
        return finish(
            Monotonicity.NOT_MONOTONE, "exact-ordering-witness", (witness, part)
        )
    return finish(
        direct, "derivative-sign", (part,) if (return_result or intervals) else ()
    )


def _classify_monotonicity(
    expr,
    var,
    *,
    domain=True,
    intervals=False,
    return_result=False,
    partition=None,
    _domain_ready=False,
):
    """Classify global monotonicity on the natural/requested real domain.

    Expensive critical-point decomposition is deliberately delayed until cheap
    registry, closure, and composition proofs have failed.  When needed, the
    shared :func:`monotonic_partition` certificate is reused by range,
    injectivity, inverse, oscillation, and extrema consumers.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if not expr.has(var):
        return Monotonicity.CONSTANT
    natural = _domain(expr, var)
    assumptions = (
        sp.sympify(domain)
        if _domain_ready
        else combine_domain_condition(natural, domain, var)
    )
    with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
        assumptions = sp.simplify(assumptions)
    if not _domain_ready:
        split = evaluate_parameter_strata(
            assumptions,
            var,
            lambda branch_domain: _classify_monotonicity(
                expr,
                var,
                domain=branch_domain,
                intervals=intervals,
                return_result=return_result,
                partition=None,
                _domain_ready=True,
            ),
            source="monotonicity-parameter-stratification",
        )
        if split is not None:
            return split
    ds = condition_to_set(assumptions, var)
    geom_domain = ds if isinstance(ds, sp.Set) else S.Reals

    early = _monotonicity_theorem_stage(
        expr,
        var,
        assumptions=assumptions,
        domain=domain,
        natural=natural,
        domain_set=ds,
        geometry_domain=geom_domain,
        intervals=intervals,
        return_result=return_result,
        partition=partition,
    )
    if early is not NO_STAGE_RESULT:
        return early

    return _monotonicity_fallback_stage(
        expr,
        var,
        assumptions=assumptions,
        geometry_domain=geom_domain,
        intervals=intervals,
        return_result=return_result,
        partition=partition,
    )


def _weaken_monotonicity(value):
    """Map strict monotonicity classes to their weak counterparts."""
    if value is Monotonicity.STRICTLY_INCREASING:
        return Monotonicity.NONDECREASING
    if value is Monotonicity.STRICTLY_DECREASING:
        return Monotonicity.NONINCREASING
    return value


def _monotonicity_implications(value):
    """Return ordinary monotonicity properties implied by ``value``."""
    if value is Monotonicity.STRICTLY_INCREASING:
        return (Monotonicity.NONDECREASING,)
    if value is Monotonicity.STRICTLY_DECREASING:
        return (Monotonicity.NONINCREASING,)
    if value is Monotonicity.CONSTANT:
        return (Monotonicity.NONDECREASING, Monotonicity.NONINCREASING)
    return ()


def _format_monotonicity(result, *, strict_inequalities, all_properties):
    """Apply public strictness and implication options to a classification."""
    weaken = (lambda value: value) if strict_inequalities else _weaken_monotonicity
    if isinstance(result, PropertyResult):
        primary = weaken(result.primary) if result.primary is not None else None
        cases = tuple(
            ConditionalExpression(weaken(case.value), case.condition)
            for case in result.cases
        )
        implied = (
            _monotonicity_implications(primary)
            if all_properties
            else result.all_properties
        )
        return PropertyResult(
            primary,
            cases=cases,
            all_properties=implied,
            source=result.source,
            exceptions=result.exceptions,
            evidence=result.evidence,
        )
    primary = weaken(result)
    if all_properties:
        return PropertyResult(
            primary,
            all_properties=_monotonicity_implications(primary),
            source="monotonicity",
        )
    return primary


def monotonicity(
    expr,
    var,
    *,
    domain=True,
    strict_inequalities=True,
    all_properties=False,
    intervals=False,
    return_result=False,
    partition=None,
):
    """Classify global monotonicity on the natural/requested real domain.

    ``strict_inequalities=False`` weakens strict increasing/decreasing results.
    ``all_properties=True`` returns the primary class together with its ordinary
    monotonicity implications. Detail and interval requests retain certificates.
    """
    result = _classify_monotonicity(
        expr,
        var,
        domain=domain,
        intervals=intervals,
        return_result=return_result,
        partition=partition,
    )
    return _format_monotonicity(
        result, strict_inequalities=strict_inequalities, all_properties=all_properties
    )
