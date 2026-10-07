from __future__ import annotations

import itertools

import sympy as sp
from sympy import S

from .. import _semialg
from ..engine import evaluate_property
from ..exactness import certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import ConditionalExpression, PropertyResult, Sign, TruthValue
from ..stratification import (
    evaluate_parameter_strata,
    implication_guard,
    normalize_cases,
    semantic_condition,
)
from ..transcendental.algebraize import algebraize
from ..transcendental.roots import transcendental_roots
from .common import (
    affine_real_parameter_guard,
    combine_domain_condition,
    condition_to_set,
    exact_order_real_points,
)
from .domain import domain as _domain
from .piecewise import combine_piecewise_sign, piecewise_branches
from .special_functions import positive_axis_sign_strata


def _constant_sign(expr):
    if expr == 0:
        return Sign.ZERO
    if expr.is_positive is True:
        return Sign.POSITIVE
    if expr.is_negative is True:
        return Sign.NEGATIVE
    if expr.is_nonnegative is True:
        return Sign.NONNEGATIVE
    if expr.is_nonpositive is True:
        return Sign.NONPOSITIVE
    return Sign.UNKNOWN


def _registry_sign(expr, var, assumptions=S.true):
    """Use catalogue sign loci as global sign proofs/parameter strata.

    Registry ``Positive`` etc. values are predicates in the function variable,
    not booleans saying that the function is globally positive.  We therefore
    prove that the requested domain implies each locus and retain any residual
    parameter-only guard.
    """
    cases = []
    for prop, val in (
        ("Positive", Sign.POSITIVE),
        ("Negative", Sign.NEGATIVE),
        ("NonNegative", Sign.NONNEGATIVE),
        ("NonPositive", Sign.NONPOSITIVE),
    ):
        r = evaluate_property(expr, prop, var=var, mode="real", assumptions=assumptions)
        rcases = r.cases or (
            () if r.primary is None else (ConditionalExpression(r.primary, S.true),)
        )
        for c in rcases:
            guard = implication_guard(assumptions, c.value, var)
            if guard is None:
                continue
            try:
                guard = sp.simplify_logic(
                    sp.And(semantic_condition(c.condition), guard)
                )
            except SYMBOLIC_QUERY_ERRORS:
                guard = sp.And(semantic_condition(c.condition), guard)
            if guard is not S.false:
                cases.append(ConditionalExpression(val, guard))
    if not cases:
        return Sign.UNKNOWN
    out = normalize_cases(cases, source="registry-sign-stratification")
    return out.primary if out.primary is not None else out


def _constant_sign_strata(expr):
    """Exact parameter strata for a constant whose real sign is unresolved."""
    expr = sp.sympify(expr)
    if not expr.free_symbols or expr.is_real is False:
        return Sign.UNKNOWN
    real_guard = S.true if expr.is_real is True else sp.Q.real(expr)
    return normalize_cases(
        (
            ConditionalExpression(Sign.POSITIVE, sp.And(real_guard, expr > 0)),
            ConditionalExpression(Sign.ZERO, sp.And(real_guard, sp.Eq(expr, 0))),
            ConditionalExpression(Sign.NEGATIVE, sp.And(real_guard, expr < 0)),
        ),
        source="constant-parameter-sign",
    )


def _semialg_sign(expr, var, assumptions):
    if not _semialg.eligible_problem(expr, var, domain=assumptions):
        return Sign.UNKNOWN
    try:
        raw = _semialg.sign(expr, var, assumptions=assumptions, return_result=True)
        value = _semialg.normalize_sign(raw)
        if value is not Sign.UNKNOWN:
            return value
        pieces = _semialg.sign_partition(expr, var, domain=assumptions)
        signs = {
            str(
                getattr(
                    piece,
                    "classification",
                    piece[0] if isinstance(piece, tuple) else piece,
                )
            )
            for piece in pieces
        }
        if "positive" in signs and "negative" in signs:
            return Sign.SIGN_CHANGING
        if signs == {"zero"}:
            return Sign.ZERO
        return Sign.UNKNOWN
    except (_semialg.SemialgUnavailable, *SYMBOLIC_QUERY_ERRORS):
        return Sign.UNKNOWN


def _multiply_sign_values(a, b):
    if a is Sign.ZERO or b is Sign.ZERO:
        return Sign.ZERO
    strict = {Sign.POSITIVE: 1, Sign.NEGATIVE: -1}
    weak = {Sign.NONNEGATIVE: 1, Sign.NONPOSITIVE: -1}
    if a in strict and b in strict:
        return Sign.POSITIVE if strict[a] * strict[b] > 0 else Sign.NEGATIVE
    if a in strict or a in weak:
        sa = strict.get(a, weak.get(a))
        if b in strict or b in weak:
            sb = strict.get(b, weak.get(b))
            return Sign.NONNEGATIVE if sa * sb > 0 else Sign.NONPOSITIVE
    return Sign.UNKNOWN


def _multiply_sign_results(values):
    cases = [ConditionalExpression(Sign.POSITIVE, S.true)]
    conditional = False
    for value in values:
        vcases = (
            value.cases
            if isinstance(value, PropertyResult)
            else (ConditionalExpression(value, S.true),)
        )
        conditional = conditional or isinstance(value, PropertyResult)
        new = []
        for a in cases:
            for b in vcases:
                v = _multiply_sign_values(a.value, b.value)
                if v is Sign.UNKNOWN:
                    continue
                try:
                    cond = sp.simplify_logic(sp.And(a.condition, b.condition))
                except SYMBOLIC_QUERY_ERRORS:
                    cond = sp.And(a.condition, b.condition)
                if cond is not S.false:
                    new.append(ConditionalExpression(v, cond))
        if not new:
            return Sign.UNKNOWN
        cases = new
    out = normalize_cases(cases, source="product-sign-stratification")
    if not conditional and out.primary is not None:
        return out.primary
    return out.primary if out.primary is not None else out


def _is_positive_axis_domain(domain, var):
    ds = condition_to_set(domain, var)
    if (
        isinstance(ds, sp.Interval)
        and ds.start == 0
        and ds.left_open
        and ds.end is S.Infinity
    ):
        return True
    try:
        return (
            sp.simplify_logic(sp.Implies(semantic_condition(domain), var > 0)) is S.true
        )
    except SYMBOLIC_QUERY_ERRORS:
        return False


def _special_sign_strata(expr, var, domain):
    if not _is_positive_axis_domain(domain, var):
        return Sign.UNKNOWN
    facts = positive_axis_sign_strata(expr, var)
    if not facts:
        return Sign.UNKNOWN
    out = normalize_cases(
        (ConditionalExpression(v, g) for v, g, *_ in facts),
        source="special-function-sign-stratification",
    )
    return out.primary if out.primary is not None else out


def _structural_sign(expr, var, domain):
    if expr.func is sp.exp:
        affine = affine_real_parameter_guard(expr.args[0], var)
        if affine is not None:
            _a, _b, guard = affine
            if guard is S.true:
                return Sign.POSITIVE
            if guard is not S.false:
                return normalize_cases(
                    (ConditionalExpression(Sign.POSITIVE, guard),),
                    source="real-exponential-positivity",
                )
        elif expr.args[0].is_real is True:
            return Sign.POSITIVE
        return Sign.UNKNOWN
    if expr.func is sp.gamma:
        ds = condition_to_set(domain, var)
        pos = (
            isinstance(ds, sp.Interval)
            and ds.start.is_nonnegative is True
            and (ds.start != 0 or ds.left_open)
        )
        pos = (
            pos
            or domain == (var > 0)
            or (isinstance(domain, sp.And) and (var > 0) in domain.args)
        )
        if pos:
            return Sign.POSITIVE
    if (
        expr.func is sp.besselk
        and len(expr.args) == 2
        and expr.args[1] == var
        and expr.args[0].is_real is True
    ):
        ds = condition_to_set(domain, var)
        pos = (
            isinstance(ds, sp.Interval)
            and ds.start.is_nonnegative is True
            and (ds.start != 0 or ds.left_open)
        )
        pos = (
            pos
            or domain == (var > 0)
            or (isinstance(domain, sp.And) and (var > 0) in domain.args)
        )
        if pos:
            return Sign.POSITIVE
    if expr.func is sp.Abs:
        # Strictness follows only if zero can be excluded.
        return Sign.NONNEGATIVE
    if expr.func is sp.Pow:
        base, e = expr.args
        if e.is_Integer and e.is_even:
            if e.is_negative:
                return Sign.POSITIVE
            return Sign.NONNEGATIVE
        if e.is_Integer and e.is_odd:
            return _function_sign(base, var, domain)
    if expr.func is sp.Mul:
        vals = [_function_sign(a, var, domain) for a in expr.args]
        return _multiply_sign_results(vals)
    return Sign.UNKNOWN


def _sample_sign(value):
    value = sp.simplify(value)
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    if certify_zero(value) is True:
        return 0
    if value.is_number:
        try:
            n = sp.N(value, 30)
            if n.is_real:
                if n > 0:
                    return 1
                if n < 0:
                    return -1
                if n == 0:
                    return 0
        except SYMBOLIC_QUERY_ERRORS:
            pass
    return None


_exact_order_real_points = exact_order_real_points


def _interval_sign(expr, var, domain):
    """Certify sign on a connected interval from continuity and exact zeros."""
    ds = condition_to_set(domain, var)
    if not isinstance(ds, sp.Interval) or ds.measure == 0:
        return Sign.UNKNOWN
    try:
        from sympy.calculus.util import continuous_domain

        if continuous_domain(expr, var, ds) != ds:
            return Sign.UNKNOWN
    except SYMBOLIC_QUERY_ERRORS:
        return Sign.UNKNOWN
    try:
        zeros = sp.solveset(sp.Eq(expr, 0), var, domain=ds)
    except SYMBOLIC_QUERY_ERRORS:
        return Sign.UNKNOWN
    if zeros is S.EmptySet:
        if ds.start is -sp.oo and ds.end is sp.oo:
            sample = S.Zero
        elif ds.start is -sp.oo:
            sample = ds.end - 1
        elif ds.end is sp.oo:
            sample = ds.start + 1
        else:
            sample = (ds.start + ds.end) / 2
        sg = _sample_sign(expr.subs(var, sample))
        return Sign.POSITIVE if sg == 1 else Sign.NEGATIVE if sg == -1 else Sign.UNKNOWN
    if not isinstance(zeros, sp.FiniteSet):
        return Sign.UNKNOWN
    interior = [z for z in zeros if z in ds and z != ds.start and z != ds.end]
    interior = _exact_order_real_points(interior)
    if interior is None:
        return Sign.UNKNOWN
    pts = [ds.start, *list(interior), ds.end]
    signs = []
    for a, b in itertools.pairwise(pts):
        if a is -sp.oo and b is sp.oo:
            sample = S.Zero
        elif a is -sp.oo:
            sample = b - 1
        elif b is sp.oo:
            sample = a + 1
        else:
            sample = (a + b) / 2
        sg = _sample_sign(expr.subs(var, sample))
        if sg is None:
            return Sign.UNKNOWN
        if sg:
            signs.append(sg)
    if signs and all(s == 1 for s in signs):
        return Sign.NONNEGATIVE
    if signs and all(s == -1 for s in signs):
        return Sign.NONPOSITIVE
    if 1 in signs and -1 in signs:
        return Sign.SIGN_CHANGING
    return Sign.UNKNOWN


def _function_sign(expr, var, domain):
    """Run the shared exact sign-proof ladder on a requested real domain."""
    expr = sp.sympify(expr)
    if isinstance(expr, sp.Piecewise):
        branches = piecewise_branches(expr, var, domain)
        if branches:
            values = [_function_sign(b.expression, var, b.condition) for b in branches]
            combined = combine_piecewise_sign(values)
            if combined is not Sign.UNKNOWN:
                return combined
    if not expr.has(var):
        s = _constant_sign(expr)
        return _constant_sign_strata(expr) if s is Sign.UNKNOWN else s
    # Real sign is undefined for an affine exponential whose exponent is
    # provably non-real on the real analysis axis.  Do not let later
    # algebraization reinterpret exp(complex) as a positive real variable.
    if expr.func is sp.exp:
        affine = affine_real_parameter_guard(expr.args[0], var)
        if affine is not None and affine[2] is S.false:
            return Sign.UNKNOWN
    if expr.func is sp.Pow and expr.exp.is_real is False:
        return Sign.UNKNOWN
    # Specialized theorem strata can be sharper than generic product closure
    # (for example a Bessel quotient whose numerator/denominator have weaker
    # individual positivity guards).
    s = _special_sign_strata(expr, var, domain)
    if s is not Sign.UNKNOWN:
        return s
    # Structural facts are globally valid on the natural domain.
    s = _structural_sign(expr, var, domain)
    if s is not Sign.UNKNOWN:
        return s
    s = _registry_sign(expr, var, domain)
    if s is not Sign.UNKNOWN:
        return s
    s = _semialg_sign(expr, var, domain)
    if s is not Sign.UNKNOWN:
        return s

    # Exact transcendental algebraization: e.g. exp(2*x)-3*exp(x)+2 is
    # a polynomial in u=exp(x), u>0.  The transformed problem can then use
    # the semialgebraic/interval proof machinery without treating exp as an
    # opaque numeric function.
    ar = algebraize(expr, var)
    if ar.changed and ar.exact and len(ar.substitutions) == 1:
        u, _ = ar.substitutions[0]
        transformed_domain = ar.constraints
        ts = _function_sign(ar.expression, u, transformed_domain)
        if ts is not Sign.UNKNOWN:
            return ts

    # Certified root partitioning on a single interval.  Unlike a numerical
    # scan, this path is used only when the root engine proves completeness.
    ds = condition_to_set(domain, var)
    if isinstance(ds, sp.Interval):
        rr = transcendental_roots(expr, var, ds, allow_numerical=False)
        if rr.complete and rr.periodic is None:
            return _interval_sign(expr, var, domain)
    return _interval_sign(expr, var, domain)


def _function_sign_public_impl(
    expr, var=None, *, domain=True, strict_inequalities=True
):
    """Classify the sign of a real expression over its natural/requested domain."""
    expr = sp.sympify(expr)
    if var is None:
        syms = tuple(expr.free_symbols)
        if len(syms) != 1:
            return _constant_sign(expr)
        var = syms[0]
    var = sp.sympify(var)
    natural = _domain(expr, var)
    assumptions = combine_domain_condition(natural, domain, var)
    split = evaluate_parameter_strata(
        assumptions,
        var,
        lambda branch_domain: _function_sign(expr, var, branch_domain),
        source="sign-parameter-stratification",
    )
    s = split if split is not None else _function_sign(expr, var, assumptions)
    if strict_inequalities:
        return s

    # Weak query can safely weaken strict classifications, including strata.
    def weaken(v):
        if v is Sign.POSITIVE:
            return Sign.NONNEGATIVE
        if v is Sign.NEGATIVE:
            return Sign.NONPOSITIVE
        return v

    if isinstance(s, PropertyResult):
        return normalize_cases(
            (ConditionalExpression(weaken(c.value), c.condition) for c in s.cases),
            primary=weaken(s.primary) if s.primary is not None else None,
            source=s.source,
            evidence=s.evidence,
        )
    return weaken(s)


def _sign_implications(value):
    """Return weak sign properties implied by an exact sign class."""
    if value is Sign.POSITIVE:
        return (Sign.NONNEGATIVE,)
    if value is Sign.NEGATIVE:
        return (Sign.NONPOSITIVE,)
    if value is Sign.ZERO:
        return (Sign.NONNEGATIVE, Sign.NONPOSITIVE)
    return ()


def sign(
    expr, var=None, *, domain=True, strict_inequalities=True, all_properties=False
):
    """Classify sign, optionally returning weaker implied sign properties."""
    expr = sp.sympify(expr)
    if var is None:
        syms = tuple(expr.free_symbols)
        if len(syms) != 1:
            result = _constant_sign(expr)
            return (
                PropertyResult(
                    result,
                    all_properties=_sign_implications(result),
                    source="constant-sign",
                )
                if all_properties
                else result
            )
        var = syms[0]
    from ..proof_cache import property_proof_context

    result = property_proof_context(expr, var, domain=domain).sign(
        strict_inequalities=strict_inequalities
    )
    if not all_properties:
        return result
    if isinstance(result, PropertyResult):
        return PropertyResult(
            result.primary,
            cases=result.cases,
            all_properties=_sign_implications(result.primary),
            source=result.source,
            exceptions=result.exceptions,
            evidence=result.evidence,
        )
    return PropertyResult(
        result, all_properties=_sign_implications(result), source="sign"
    )


def prove_relation(relation, var=None, *, domain=True):
    """Prove a scalar real relation by reducing it to the shared sign engine.

    Equality and order relations are normalized to the sign of ``lhs-rhs``.
    Conditional parameter results are preserved
    unknown Boolean.
    """
    relation = sp.sympify(relation)
    if relation in (S.true, S.false):
        return TruthValue.TRUE if relation is S.true else TruthValue.FALSE
    if not isinstance(relation, sp.core.relational.Relational):
        raise TypeError("relation must be a SymPy equality or inequality")
    expr = sp.simplify(relation.lhs - relation.rhs)
    sign_value = sign(expr, var, domain=domain)

    def classify(value):
        if isinstance(relation, sp.StrictGreaterThan):
            return (
                TruthValue.TRUE
                if value is Sign.POSITIVE
                else (
                    TruthValue.FALSE
                    if value in {Sign.ZERO, Sign.NONPOSITIVE, Sign.NEGATIVE}
                    else TruthValue.UNKNOWN
                )
            )
        if isinstance(relation, sp.GreaterThan):
            return (
                TruthValue.TRUE
                if value in {Sign.POSITIVE, Sign.NONNEGATIVE, Sign.ZERO}
                else (
                    TruthValue.FALSE if value is Sign.NEGATIVE else TruthValue.UNKNOWN
                )
            )
        if isinstance(relation, sp.StrictLessThan):
            return (
                TruthValue.TRUE
                if value is Sign.NEGATIVE
                else (
                    TruthValue.FALSE
                    if value in {Sign.ZERO, Sign.NONNEGATIVE, Sign.POSITIVE}
                    else TruthValue.UNKNOWN
                )
            )
        if isinstance(relation, sp.LessThan):
            return (
                TruthValue.TRUE
                if value in {Sign.NEGATIVE, Sign.NONPOSITIVE, Sign.ZERO}
                else (
                    TruthValue.FALSE if value is Sign.POSITIVE else TruthValue.UNKNOWN
                )
            )
        if isinstance(relation, sp.Equality):
            return (
                TruthValue.TRUE
                if value is Sign.ZERO
                else (
                    TruthValue.FALSE
                    if value in {Sign.POSITIVE, Sign.NEGATIVE, Sign.SIGN_CHANGING}
                    else TruthValue.UNKNOWN
                )
            )
        if isinstance(relation, sp.Unequality):
            return (
                TruthValue.FALSE
                if value is Sign.ZERO
                else (
                    TruthValue.TRUE
                    if value in {Sign.POSITIVE, Sign.NEGATIVE}
                    else TruthValue.UNKNOWN
                )
            )
        return TruthValue.UNKNOWN

    if isinstance(sign, PropertyResult):
        return normalize_cases(
            (ConditionalExpression(classify(c.value), c.condition) for c in sign.cases),
            primary=classify(sign.primary) if sign.primary is not None else None,
            source="relation-from-sign",
            evidence=sign.evidence,
        )
    return classify(sign_value)
