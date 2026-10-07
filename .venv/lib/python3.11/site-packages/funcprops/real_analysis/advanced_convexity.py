from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import ConditionalExpression, Convexity, Monotonicity, Sign
from ..stratification import implication_guard, normalize_cases
from .common import affine_real_parameter_guard, condition_to_set
from .monotonicity import monotonicity
from .sign import sign


@dataclass(frozen=True)
class ConvexityCertificate:
    """Replayable evidence for one ordinary convexity classification."""

    classification: Convexity
    domain: object
    method: str
    modulus: sp.Expr | None = None
    proof: tuple = ()


def parameterized_power_convexity(expr, var, domain):
    """Exact parameter strata for x**a on a positive real domain."""
    if expr.func is not sp.Pow or expr.base != var:
        return None
    a = expr.exp
    if not (a.free_symbols - {var}) and a.is_number:
        return None
    ds = condition_to_set(domain, var)
    # This stratification is valid only where x>0.
    positive_domain = False
    if isinstance(ds, sp.Interval):
        positive_domain = (ds.start.is_nonnegative is True and ds.start != 0) or (
            ds.start == 0 and ds.left_open
        )
    if not positive_domain:
        # Rich parameter domains can contain Contains/NotContains predicates,
        # which SymPy's LRA-backed ask() does not accept.  Use the package's
        # three-valued implication reducer
        # formulae into ask().
        positive_domain = implication_guard(domain, var > 0, var) is S.true
    if not positive_domain and domain != (var > 0):
        return None
    real_a = S.true if a.is_real is True else sp.Q.real(a)
    cases = (
        ConditionalExpression(
            Convexity.STRICTLY_CONVEX, sp.And(real_a, sp.Or(a < 0, a > 1))
        ),
        ConditionalExpression(Convexity.STRICTLY_CONCAVE, sp.And(real_a, a > 0, a < 1)),
        ConditionalExpression(Convexity.AFFINE, sp.Or(sp.Eq(a, 0), sp.Eq(a, 1))),
    )
    return normalize_cases(cases, source="power-parameter-stratification")


def parameterized_exp_convexity(expr, var):
    if expr.func is not sp.exp:
        return None
    affine = affine_real_parameter_guard(expr.args[0], var)
    if affine is None:
        return None
    a, _b, real_guard = affine
    if real_guard is S.false:
        return None
    if not a.free_symbols and real_guard is S.true:
        return None
    return normalize_cases(
        (
            ConditionalExpression(
                Convexity.STRICTLY_CONVEX, sp.And(real_guard, sp.Ne(a, 0))
            ),
            ConditionalExpression(Convexity.AFFINE, sp.And(real_guard, sp.Eq(a, 0))),
        ),
        source="exp-affine-parameter-stratification",
    )


def strong_convexity(expr, var, domain):
    """Prove strong convexity/concavity by an exact second-derivative bound."""
    try:
        from ..proof_cache import cached_derivative

        d2 = sp.simplify(cached_derivative(expr, var, order=2))
    except SYMBOLIC_QUERY_ERRORS:
        return None
    if d2 == 0:
        return None
    ds = domain if isinstance(domain, sp.Set) else condition_to_set(domain, var)
    if not isinstance(ds, sp.Interval):
        return None
    # Constant Hessian is the clean exact case.
    if not d2.has(var):
        if d2.is_positive is True:
            return Convexity.STRONGLY_CONVEX, ConvexityCertificate(
                Convexity.STRONGLY_CONVEX, ds, "second-derivative-bound", d2, (d2,)
            )
        if d2.is_negative is True:
            return Convexity.STRONGLY_CONCAVE, ConvexityCertificate(
                Convexity.STRONGLY_CONCAVE, ds, "second-derivative-bound", -d2, (d2,)
            )
    # Exact infimum/supremum on an interval. Avoid numerical optimization.
    try:
        from sympy.calculus.util import maximum, minimum

        lo = minimum(d2, var, ds)
        if lo is not None and lo is not -sp.oo and sp.ask(sp.Q.positive(lo)) is True:
            return Convexity.STRONGLY_CONVEX, ConvexityCertificate(
                Convexity.STRONGLY_CONVEX,
                ds,
                "second-derivative-bound",
                sp.simplify(lo),
                (d2, lo),
            )
        hi = maximum(d2, var, ds)
        if hi is not None and hi is not sp.oo and sp.ask(sp.Q.negative(hi)) is True:
            return Convexity.STRONGLY_CONCAVE, ConvexityCertificate(
                Convexity.STRONGLY_CONCAVE,
                ds,
                "second-derivative-bound",
                sp.simplify(-hi),
                (d2, hi),
            )
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return None


def convexity_implications(primary):
    """Return ordinary convexity/concavity properties implied by ``primary``.

    The ordering lives in one place so default strongest-class results and
    ``all_properties=True`` cannot drift apart across scalar and multivariate
    code paths.
    """
    table = {
        Convexity.STRONGLY_CONVEX: (Convexity.STRICTLY_CONVEX, Convexity.CONVEX),
        Convexity.STRICTLY_CONVEX: (Convexity.CONVEX,),
        Convexity.CONVEX: (),
        Convexity.STRONGLY_CONCAVE: (Convexity.STRICTLY_CONCAVE, Convexity.CONCAVE),
        Convexity.STRICTLY_CONCAVE: (Convexity.CONCAVE,),
        Convexity.CONCAVE: (),
        Convexity.AFFINE: (Convexity.CONVEX, Convexity.CONCAVE),
    }
    return table.get(primary, ())


def secondary_convexity_properties(expr, var, domain, primary):
    props = list(convexity_implications(primary))
    evidence = []
    # Positive f: log-convexity/concavity reduces exactly to convexity of log(f).
    s = sign(expr, var, domain=domain)
    if s is Sign.POSITIVE:
        try:
            from ..proof_cache import cached_derivative

            d2log = sp.simplify(cached_derivative(sp.log(expr), var, order=2))
            ls = sign(d2log, var, domain=domain)
            if ls in (Sign.POSITIVE, Sign.NONNEGATIVE, Sign.ZERO):
                props.append(Convexity.LOG_CONVEX)
            if ls in (Sign.NEGATIVE, Sign.NONPOSITIVE, Sign.ZERO):
                props.append(Convexity.LOG_CONCAVE)
        except SYMBOLIC_QUERY_ERRORS:
            pass
    # Every monotone univariate function is both quasiconvex and quasiconcave.
    m = monotonicity(expr, var, domain=domain)
    if m in (
        Monotonicity.STRICTLY_INCREASING,
        Monotonicity.NONDECREASING,
        Monotonicity.STRICTLY_DECREASING,
        Monotonicity.NONINCREASING,
        Monotonicity.CONSTANT,
    ):
        props.extend((Convexity.QUASICONVEX, Convexity.QUASICONCAVE))
    # A differentiable function with derivative of one strict sign is both
    # pseudoconvex and pseudoconcave on an interval.
    try:
        ds = condition_to_set(domain, var)
        if isinstance(ds, sp.Interval):
            from ..proof_cache import cached_derivative

            sg = sign(cached_derivative(expr, var, order=1), var, domain=domain)
            if sg in (Sign.POSITIVE, Sign.NEGATIVE):
                props.extend((Convexity.PSEUDOCONVEX, Convexity.PSEUDOCONCAVE))
    except SYMBOLIC_QUERY_ERRORS:
        pass
    # deterministic order / dedupe
    out = []
    for p in props:
        if p not in out and p is not primary:
            out.append(p)
    return tuple(out), tuple(evidence)


def _positive_domain(domain, var):
    ds = condition_to_set(domain, var)
    if isinstance(ds, sp.Interval):
        return ds.start.is_nonnegative is True and (ds.start != 0 or ds.left_open)
    rel = var > 0
    return domain == rel or (isinstance(domain, sp.And) and rel in domain.args)


def special_function_convexity(expr, var, domain):
    condition_to_set(domain, var)
    # Gamma is positive and strictly log-convex on (0, infinity); positivity
    # plus strict log-convexity also implies strict ordinary convexity.
    if expr.func is sp.gamma and expr.args == (var,) and _positive_domain(domain, var):
        return (Convexity.STRICTLY_CONVEX, Convexity.CONVEX, Convexity.LOG_CONVEX)
    # Modified Bessel K_nu(x): for real fixed order, catalogue facts already
    # give decreasing/convex; expose that as a strong special-function hint.
    if expr.func is sp.besselk and len(expr.args) == 2 and expr.args[1] == var:
        nu = expr.args[0]
        if nu.is_real is True and _positive_domain(domain, var):
            return (Convexity.STRICTLY_CONVEX, Convexity.CONVEX, Convexity.LOG_CONVEX)
    return ()
