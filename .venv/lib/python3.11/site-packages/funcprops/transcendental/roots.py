from __future__ import annotations

import contextlib
import math
from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..exactness import certify_identity_equal, certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from .analytic import analytic_neighborhood
from .periodic import PeriodicSet
from .real_sign import robust_real_sign
from .zero_geometry import ZeroGeometryCertificate, zero_geometry


@dataclass(frozen=True)
class RootCertificate:
    domain: object
    root_count: int | None
    method: str
    proof: tuple[object, ...] = ()
    complete: bool = True
    multiplicity: int | None = None


@dataclass(frozen=True)
class CertifiedRootInterval:
    interval: sp.Interval
    certificate: RootCertificate


@dataclass(frozen=True)
class ParametricZeroFamily:
    count: sp.Expr
    support: sp.Set
    guard: sp.Expr = S.true
    simple: bool | None = None
    interlaces_prev_degree: bool | None = None
    interlacing: tuple[str, ...] = ()
    source: str = ""


@dataclass(frozen=True)
class RootIsolationResult:
    roots: tuple[sp.Expr, ...] = ()
    zero_intervals: tuple[sp.Interval, ...] = ()
    isolated_intervals: tuple[CertifiedRootInterval, ...] = ()
    periodic: PeriodicSet | None = None
    complete: bool = False
    method: str = "unknown"
    evidence: tuple[object, ...] = ()
    certificates: tuple[RootCertificate, ...] = ()
    parametric_families: tuple[ParametricZeroFamily, ...] = ()
    geometry: tuple[ZeroGeometryCertificate, ...] = ()

    @property
    def certified(self):
        return self.complete


def _sort_unique(vals):
    out = []
    for v in sorted((sp.simplify(q) for q in vals), key=sp.default_sort_key):
        if not any(certify_identity_equal(v, u) is True for u in out):
            out.append(v)
    return tuple(out)


def _affine_data(arg, var):
    try:
        p = sp.Poly(sp.expand(arg), var)
        if p.degree() != 1:
            return None
        return p.nth(1), p.nth(0)
    except SYMBOLIC_QUERY_ERRORS:
        return None


def _periodic(expr, var, domain):
    if expr.func not in (sp.sin, sp.cos, sp.tan):
        return None
    ab = _affine_data(expr.args[0], var)
    if ab is None:
        return None
    a, b = ab
    if a == 0 or a.is_real is False:
        return None
    target = S.Zero if expr.func in (sp.sin, sp.tan) else sp.pi / 2
    p = sp.simplify(sp.pi / sp.Abs(a))
    if p.is_positive is not True:
        return None
    x0 = sp.simplify((target - b) / a)
    base = sp.FiniteSet(sp.simplify(sp.Mod(x0, p)))
    ps = PeriodicSet(base, p, var)
    if domain is S.Reals:
        c = RootCertificate(
            S.Reals, None, "periodic_formula", (sp.Eq(expr, 0), ps), True
        )
        return RootIsolationResult(
            periodic=ps, complete=True, method="structural_periodic", certificates=(c,)
        )
    return None


def _special(expr, var, domain):
    if expr.func in (sp.erf, sp.erfi) and len(expr.args) == 1:
        ab = _affine_data(expr.args[0], var)
        if ab:
            a, b = ab
            r = sp.simplify(-b / a)
            rts = () if sp.Contains(r, domain) is S.false else (r,)
            c = RootCertificate(
                domain,
                len(rts),
                "special_function_zero",
                ("erf/erfi has unique real zero at 0",),
                True,
            )
            return RootIsolationResult(
                rts, complete=True, method="special", certificates=(c,)
            )
    if expr.func is sp.gamma:
        c = RootCertificate(
            domain, 0, "special_function_no_zeros", ("gamma has no zeros",), True
        )
        return RootIsolationResult(
            (), complete=True, method="special_no_zeros", certificates=(c,)
        )
    return None


def _guard_proved(g):
    g = sp.sympify(g)
    if g is S.true:
        return True
    try:
        q = sp.refine(g)
        if q is S.true:
            return True
    except SYMBOLIC_QUERY_ERRORS:
        pass
    try:
        q = sp.simplify(g)
        if q is S.true:
            return True
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return False


def _orthogonal_special(expr, var, domain):
    """Return exact zero-family metadata for audited orthogonal polynomials.

    Parameter guards are retained symbolically; completeness is asserted only
    when the support and degree hypotheses are proved on the requested domain.
    """
    # Formula-level DLMF zero distributions.  For unresolved parameter guards
    # retain exact conditional metadata without claiming completeness.
    if not expr.args or expr.args[-1] != var:
        return None
    f = expr.func
    n = expr.args[0] if expr.args else None
    n_guard = sp.And(sp.Contains(n, S.Integers), n >= 1) if n is not None else S.false
    support = None
    guard = n_guard
    source = ""
    interlace = True
    if f is sp.jacobi and len(expr.args) == 4:
        a, b = expr.args[1:3]
        support = sp.Interval.open(-1, 1)
        guard = sp.And(n_guard, a > -1, b > -1)
        source = "DLMF 18.16.E1 + 18.2(vi)"
    elif f is sp.gegenbauer and len(expr.args) == 3:
        lam = expr.args[1]
        support = sp.Interval.open(-1, 1)
        guard = sp.And(n_guard, lam > -sp.Rational(1, 2), sp.Ne(lam, 0))
        source = "DLMF 18.16(iii) + 18.2(vi)"
    elif f is sp.assoc_laguerre and len(expr.args) == 3:
        a = expr.args[1]
        support = sp.Interval.open(0, sp.oo)
        guard = sp.And(n_guard, a > -1)
        source = "DLMF 18.16.E9 + 18.2(vi)"
    elif f is sp.laguerre and len(expr.args) == 2:
        support = sp.Interval.open(0, sp.oo)
        source = "DLMF 18.16.E9 + 18.2(vi)"
    elif f in tuple(
        q
        for q in (getattr(sp, "hermite", None), getattr(sp, "hermite_prob", None))
        if q is not None
    ):
        support = S.Reals
        source = "DLMF 18.16(v) + 18.2(vi)"
    elif f is sp.chebyshevt and len(expr.args) == 2:
        support = sp.Interval.open(-1, 1)
        source = "DLMF 18.5.E1 + 18.16(iii)"
        if _guard_proved(n_guard):
            k = sp.Dummy("k", integer=True)
            roots = sp.ImageSet(
                sp.Lambda(k, sp.cos((2 * k - 1) * sp.pi / (2 * n))), sp.Range(1, n + 1)
            )
            with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
                roots = sp.Intersection(roots, domain)
            c = RootCertificate(domain, None, "DLMF_chebyshev_formula", (source,), True)
            fam = ParametricZeroFamily(n, support, n_guard, True, True, (), source)
            return RootIsolationResult(
                complete=True,
                method="dlmf_indexed_zeros",
                evidence=(roots,),
                certificates=(c,),
                parametric_families=(fam,),
            )
    elif f is sp.chebyshevu and len(expr.args) == 2:
        support = sp.Interval.open(-1, 1)
        source = "DLMF 18.5.E2 + 18.16(iii)"
        if _guard_proved(n_guard):
            k = sp.Dummy("k", integer=True)
            roots = sp.ImageSet(
                sp.Lambda(k, sp.cos(k * sp.pi / (n + 1))), sp.Range(1, n + 1)
            )
            with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
                roots = sp.Intersection(roots, domain)
            c = RootCertificate(domain, None, "DLMF_chebyshev_formula", (source,), True)
            fam = ParametricZeroFamily(n, support, n_guard, True, True, (), source)
            return RootIsolationResult(
                complete=True,
                method="dlmf_indexed_zeros",
                evidence=(roots,),
                certificates=(c,),
                parametric_families=(fam,),
            )
    elif f is sp.assoc_legendre and len(expr.args) == 3:
        mu = expr.args[1]
        support = sp.Interval.open(-1, 1)
        guard = sp.And(
            sp.Contains(n, S.Naturals0), sp.Contains(mu, S.Integers), sp.Abs(mu) <= n
        )
        count = sp.simplify(n - sp.Abs(mu))
        source = "DLMF 14.16(ii), integer-degree case"
        fam = ParametricZeroFamily(count, support, guard, None, None, (), source)
        contained = domain == support
        if _guard_proved(guard) and contained:
            c = RootCertificate(
                domain,
                None,
                "DLMF_Ferrers_zero_count",
                (source, ("count", count), ("support", support)),
                True,
            )
            return RootIsolationResult(
                complete=True,
                method="dlmf_zero_distribution",
                evidence=(sp.ConditionSet(var, sp.Eq(expr, 0), support),),
                certificates=(c,),
                parametric_families=(fam,),
            )
        return RootIsolationResult(
            complete=False,
            method="conditional_dlmf_zero_distribution",
            evidence=(sp.ConditionSet(var, sp.Eq(expr, 0), support),),
            parametric_families=(fam,),
        )
    elif f in (sp.jn, sp.yn) and len(expr.args) == 2:
        support = sp.Interval.open(0, sp.oo)
        guard = sp.And(sp.Contains(n, S.Integers), n >= 0)
        source = "DLMF 10.58.E1 + 10.21.E2"
        rel = ("positive zeros interlace with adjacent order",)
        fam = ParametricZeroFamily(sp.oo, support, guard, True, True, rel, source)
        contained = domain == support
        if _guard_proved(guard) and contained:
            c = RootCertificate(
                domain,
                None,
                "DLMF_spherical_bessel_zeros",
                (source, ("count", sp.oo), ("support", support)),
                True,
            )
            return RootIsolationResult(
                complete=True,
                method="dlmf_zero_distribution",
                evidence=(sp.ConditionSet(var, sp.Eq(expr, 0), support),),
                certificates=(c,),
                parametric_families=(fam,),
            )
        return RootIsolationResult(
            complete=False,
            method="conditional_dlmf_zero_distribution",
            evidence=(sp.ConditionSet(var, sp.Eq(expr, 0), support),),
            parametric_families=(fam,),
        )
    if support is None:
        return None
    fam = ParametricZeroFamily(n, support, guard, True, interlace, (), source)
    contained = domain is S.Reals or domain == support
    if _guard_proved(guard) and contained:
        c = RootCertificate(
            domain,
            None,
            "DLMF_zero_distribution",
            (source, ("count", n), ("support", support)),
            True,
        )
        return RootIsolationResult(
            complete=True,
            method="dlmf_zero_distribution",
            evidence=(sp.ConditionSet(var, sp.Eq(expr, 0), support),),
            certificates=(c,),
            parametric_families=(fam,),
        )
    return RootIsolationResult(
        complete=False,
        method="conditional_dlmf_zero_distribution",
        evidence=(sp.ConditionSet(var, sp.Eq(expr, 0), support),),
        parametric_families=(fam,),
    )


def _concrete_zero_set_parts(solution):
    """Split a concrete real solution set into exact points and zero intervals.

    ``None`` means the set still contains symbolic/conditional structure and cannot
    support a completeness claim.
    """
    if solution is S.EmptySet:
        return (), ()
    if isinstance(solution, sp.FiniteSet):
        return _sort_unique(solution), ()
    if isinstance(solution, sp.Interval):
        return (), (solution,)
    if isinstance(solution, sp.Union):
        points = []
        intervals = []
        for part in solution.args:
            split = _concrete_zero_set_parts(part)
            if split is None:
                return None
            p, i = split
            points.extend(p)
            intervals.extend(i)
        return _sort_unique(points), tuple(
            sorted(set(intervals), key=sp.default_sort_key)
        )
    return None


def _identically_zero_interval(expr, var, domain):
    """Prove that a complete interval cell belongs to the zero set."""
    if not isinstance(domain, sp.Interval):
        return None
    try:
        reduced = sp.simplify(expr)
    except SYMBOLIC_QUERY_ERRORS:
        reduced = expr
    if not isinstance(reduced, sp.Piecewise) and certify_zero(reduced) is True:
        cert = RootCertificate(
            domain, None, "identically_zero_interval", (("simplified", reduced),), True
        )
        return RootIsolationResult(
            zero_intervals=(domain,),
            complete=True,
            method="identically_zero_interval",
            certificates=(cert,),
        )
    if isinstance(expr, sp.Piecewise):
        try:
            from ..real_analysis.common import condition_to_set

            exclusive = sp.piecewise_exclusive(expr, deep=True)
            for branch, condition in exclusive.args:
                branch_set = condition_to_set(condition, var)
                if branch_set is None:
                    continue
                if domain.is_subset(branch_set) is True:
                    value = sp.simplify(branch)
                    if certify_zero(value) is True:
                        cert = RootCertificate(
                            domain,
                            None,
                            "solveset_identically_zero_interval",
                            (("branch_condition", condition),),
                            True,
                        )
                        return RootIsolationResult(
                            zero_intervals=(domain,),
                            complete=True,
                            method="identically_zero_interval",
                            certificates=(cert,),
                        )
        except SYMBOLIC_QUERY_ERRORS:
            pass
    return None


def _solveset(expr, var, domain):
    try:
        sol = sp.solveset(sp.Eq(expr, 0), var, domain=domain)
    except SYMBOLIC_QUERY_ERRORS:
        return None
    concrete = _concrete_zero_set_parts(sol)
    if concrete is None:
        return None
    roots, intervals = concrete
    certs = []
    if intervals:
        certs.extend(
            RootCertificate(
                interval,
                None,
                "solveset_identically_zero_interval",
                (("solution_set", sol),),
                True,
            )
            for interval in intervals
        )
    if roots:
        certs.append(
            RootCertificate(
                sp.FiniteSet(*roots),
                len(roots),
                "solveset",
                (sol,),
                True,
            )
        )
    if not roots and not intervals:
        certs.append(RootCertificate(domain, 0, "solveset", (sol,), True))
    method = "solveset_zero_geometry" if intervals else "solveset"
    return RootIsolationResult(
        roots,
        intervals,
        complete=True,
        method=method,
        certificates=tuple(certs),
    )


def _endpoint_value(expr, var, x, side):
    from ..proof_cache import cached_endpoint_limit

    return cached_endpoint_limit(expr, var, x, direction=side, attained=False)


def _strict_derivative_sign(d, var, interval):
    # exact positivity/negativity via solveset of the negating condition
    try:
        bad = sp.solveset(d <= 0, var, domain=interval)
        if bad is S.EmptySet:
            return 1
        bad = sp.solveset(d >= 0, var, domain=interval)
        if bad is S.EmptySet:
            return -1
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return 0


def _certified_monotone_interval(expr, var, domain):
    """Certify at most one root on an interval from a strict derivative sign."""
    if not isinstance(domain, sp.Interval):
        return None
    try:
        from ..proof_cache import cached_derivative

        d = cached_derivative(expr, var, order=1)
    except SYMBOLIC_QUERY_ERRORS:
        return None
    s = _strict_derivative_sign(d, var, domain)
    if not s:
        return None
    ana = (
        analytic_neighborhood(expr, var, domain)
        if domain.start not in (-sp.oo, sp.oo) and domain.end not in (-sp.oo, sp.oo)
        else None
    )
    # real differentiability + strict derivative sign suffices; complex cert is extra evidence
    left = _endpoint_value(expr, var, domain.start, "+")
    right = _endpoint_value(expr, var, domain.end, "-")
    if left is None or right is None:
        return None
    sl = sp.sign(left)
    sr = sp.sign(right)
    cert = RootCertificate(
        domain,
        None,
        "strict_monotonicity",
        (("derivative_sign", s), ("endpoint_values", left, right), ana),
        True,
    )
    if certify_zero(left) is True and not domain.left_open:
        cert = RootCertificate(domain, 1, "strict_monotonicity", cert.proof, True)
        return RootIsolationResult(
            (domain.start,),
            complete=True,
            method="certified_monotonic",
            certificates=(cert,),
        )
    if certify_zero(right) is True and not domain.right_open:
        cert = RootCertificate(domain, 1, "strict_monotonicity", cert.proof, True)
        return RootIsolationResult(
            (domain.end,),
            complete=True,
            method="certified_monotonic",
            certificates=(cert,),
        )
    if (sl.is_positive and sr.is_positive) or (sl.is_negative and sr.is_negative):
        cert = RootCertificate(domain, 0, "strict_monotonicity", cert.proof, True)
        return RootIsolationResult(
            (), complete=True, method="certified_monotonic", certificates=(cert,)
        )
    if (sl.is_negative and sr.is_positive) or (sl.is_positive and sr.is_negative):
        # Unique root is certified even if there is no closed form. Produce a nested rational isolating interval.
        if domain.start in (-sp.oo, sp.oo) or domain.end in (-sp.oo, sp.oo):
            return None
        a = sp.Rational(domain.start) if domain.start.is_Rational else domain.start
        b = sp.Rational(domain.end) if domain.end.is_Rational else domain.end
        # exact bisection using symbolic signs; stop at 24 certified subdivisions or earlier unknown
        aa, bb = a, b
        sa = sl
        for _ in range(24):
            m = sp.simplify((aa + bb) / 2)
            sm = sp.sign(sp.simplify(expr.subs(var, m)))
            if certify_zero(sp.simplify(expr.subs(var, m))) is True:
                cert = RootCertificate(
                    domain,
                    1,
                    "strict_monotonicity+exact_bisection",
                    (*cert.proof, ("exact_root", m)),
                    True,
                )
                return RootIsolationResult(
                    (m,),
                    complete=True,
                    method="certified_monotonic",
                    certificates=(cert,),
                )
            if sm.is_positive is not True and sm.is_negative is not True:
                break
            if (sa.is_negative and sm.is_positive) or (
                sa.is_positive and sm.is_negative
            ):
                bb = m
            else:
                aa = m
                sa = sm
        iv = sp.Interval(aa, bb)
        cert = RootCertificate(
            domain,
            1,
            "strict_monotonicity+isolation",
            (*cert.proof, ("isolating_interval", iv)),
            True,
        )
        return RootIsolationResult(
            isolated_intervals=(CertifiedRootInterval(iv, cert),),
            complete=True,
            method="certified_monotonic_interval",
            certificates=(cert,),
        )
    return None


def _exact_sign(value):
    """Return -1, 0, or 1 only when the real sign is certified."""
    return robust_real_sign(value)


def _certified_multiplicity(expr, var, root, *, max_order=8):
    """Return a proved finite multiplicity, or ``None`` if not established."""
    current = expr
    for order in range(max_order + 1):
        try:
            value = sp.simplify(current.subs(var, root))
        except SYMBOLIC_QUERY_ERRORS:
            return None
        sign = _exact_sign(value)
        if order == 0 and sign not in (0, None):
            return None
        if order > 0 and sign in (-1, 1):
            return order
        if sign is None:
            return None
        try:
            current = sp.diff(current, var)
        except SYMBOLIC_QUERY_ERRORS:
            return None
    return None


def _dyadic_interior_points(interval, *, max_level=24):
    """Yield exact dyadic points approaching both ends of a finite interval."""
    if not isinstance(interval, sp.Interval):
        return
    try:
        left = float(sp.N(interval.start, 40))
        right = float(sp.N(interval.end, 40))
    except SYMBOLIC_QUERY_ERRORS:
        return
    if not (math.isfinite(left) and math.isfinite(right) and left < right):
        return
    seen = set()
    for level in range(4, max_level + 1):
        den = 1 << level
        lo_num = math.floor(left * den) + 1
        hi_num = math.ceil(right * den) - 1
        for num in (lo_num, hi_num):
            q = sp.Rational(num, den)
            if q in seen:
                continue
            seen.add(q)
            try:
                inside = interval.contains(q)
            except SYMBOLIC_QUERY_ERRORS:
                continue
            if inside is S.true:
                yield q


def _rational_sign_bracket(expr, var, interval, left_sign, right_sign):
    """Find rational points inside ``interval`` with the requested opposite signs."""
    left = None
    right = None
    for q in _dyadic_interior_points(interval):
        sign = _exact_sign(expr.subs(var, q))
        if sign == left_sign and (left is None or q < left):
            left = q
        if sign == right_sign and (right is None or q > right):
            right = q
        if left is not None and right is not None and left < right:
            return left, right
    return None


def _bisect_certified_bracket(expr, var, a, b, sa, *, steps=32):
    """Tighten an exact opposite-sign rational bracket without losing certification."""
    aa, bb = sp.Rational(a), sp.Rational(b)
    left_sign = sa
    for _ in range(steps):
        m = (aa + bb) / 2
        sm = _exact_sign(expr.subs(var, m))
        if sm == 0:
            return m, None
        if sm not in (-1, 1):
            break
        if sm == left_sign:
            aa = m
        else:
            bb = m
    return None, sp.Interval(aa, bb)


def _certified_range_on_interval(expr, var, interval):
    """Return an exact range set on ``interval`` when funcprops can certify it."""
    try:
        from ..real_analysis.common import condition_to_set, set_to_condition
        from ..real_analysis.range import range

        y = sp.Dummy("y", real=True)
        cond = range(
            expr,
            var,
            y,
            constraints=set_to_condition(interval, var),
        )
        if hasattr(cond, "cases") and hasattr(cond, "primary"):
            return None
        result = condition_to_set(cond, y)
    except (ImportError, *SYMBOLIC_QUERY_ERRORS):
        return None
    return result if isinstance(result, sp.Set) else None


def _range_sign(range_set):
    """Return a strict common sign for a certified range set, if available."""
    try:
        if sp.Contains(0, range_set) is not S.false:
            return 0
        positive = sp.Intersection(range_set, sp.Interval.open(0, sp.oo))
        negative = sp.Intersection(range_set, sp.Interval.open(-sp.oo, 0))
    except SYMBOLIC_QUERY_ERRORS:
        return 0
    if positive is not S.EmptySet and negative is S.EmptySet:
        return 1
    if negative is not S.EmptySet and positive is S.EmptySet:
        return -1
    return 0


def _certified_interval_resolution(expr, var, domain):
    """Resolve a finite real interval from exact range/derivative-range evidence.

    This is the funcprops analogue of an interval-arithmetic undecided-cell pass:
    an exact image excluding zero proves a cell root-free; an exact derivative image
    of one strict sign proves monotonicity and reduces the zero question to endpoint
    signs plus rational isolation.
    """
    if (
        not isinstance(domain, sp.Interval)
        or domain.start in (-sp.oo, sp.oo)
        or domain.end in (-sp.oo, sp.oo)
    ):
        return None

    try:
        from .real_interval import certified_real_interval_evaluate

        interval_value = certified_real_interval_evaluate(expr, var, domain)
    except ImportError:
        interval_value = None
    if interval_value is not None and interval_value.excludes_zero:
        cert = RootCertificate(
            domain,
            0,
            "certified_interval_range_exclusion",
            (("interval_evaluation", interval_value),),
            True,
        )
        return RootIsolationResult(
            complete=True,
            method="certified_interval_resolution",
            evidence=(interval_value,),
            certificates=(cert,),
        )

    image = _certified_range_on_interval(expr, var, domain)
    if image is not None:
        try:
            excludes_zero = sp.Contains(0, image) is S.false
        except SYMBOLIC_QUERY_ERRORS:
            excludes_zero = False
        if excludes_zero:
            cert = RootCertificate(
                domain,
                0,
                "certified_interval_range_exclusion",
                (("range", image),),
                True,
            )
            return RootIsolationResult(
                complete=True,
                method="certified_interval_resolution",
                evidence=(("range", image),),
                certificates=(cert,),
            )

    try:
        derivative = sp.diff(expr, var)
    except SYMBOLIC_QUERY_ERRORS:
        return None
    try:
        derivative_interval = certified_real_interval_evaluate(derivative, var, domain)
    except (NameError, ValueError):
        derivative_interval = None
    derivative_image = _certified_range_on_interval(derivative, var, domain)
    derivative_sign = (
        derivative_interval.strict_sign
        if derivative_interval is not None and derivative_interval.excludes_zero
        else _range_sign(derivative_image)
        if derivative_image is not None
        else 0
    )
    if not derivative_sign:
        return None

    left = _endpoint_value(expr, var, domain.start, "+")
    right = _endpoint_value(expr, var, domain.end, "-")
    if left is None or right is None:
        return None
    sl = _exact_sign(left)
    sr = _exact_sign(right)
    if sl is None or sr is None:
        return None
    proof = (
        ("derivative_interval", derivative_interval),
        ("derivative_range", derivative_image),
        ("derivative_sign", derivative_sign),
        ("endpoint_values", left, right),
    )

    if sl == 0 and not domain.left_open:
        mult = _certified_multiplicity(expr, var, domain.start)
        cert = RootCertificate(
            domain, 1, "interval_resolution_endpoint_root", proof, True, mult
        )
        return RootIsolationResult(
            (domain.start,),
            complete=True,
            method="certified_interval_resolution",
            certificates=(cert,),
        )
    if sr == 0 and not domain.right_open:
        mult = _certified_multiplicity(expr, var, domain.end)
        cert = RootCertificate(
            domain, 1, "interval_resolution_endpoint_root", proof, True, mult
        )
        return RootIsolationResult(
            (domain.end,),
            complete=True,
            method="certified_interval_resolution",
            certificates=(cert,),
        )
    if sl == sr and sl in (-1, 1):
        cert = RootCertificate(
            domain, 0, "interval_resolution_monotone_same_sign", proof, True
        )
        return RootIsolationResult(
            complete=True,
            method="certified_interval_resolution",
            certificates=(cert,),
        )
    if sl * sr != -1:
        return None

    bracket = _rational_sign_bracket(expr, var, domain, sl, sr)
    if bracket is None:
        return None
    a, b = bracket
    exact_root, isolating = _bisect_certified_bracket(expr, var, a, b, sl)
    if exact_root is not None:
        mult = _certified_multiplicity(expr, var, exact_root)
        cert = RootCertificate(
            domain,
            1,
            "interval_resolution_exact_bisection",
            (*proof, ("rational_bracket", a, b), ("exact_root", exact_root)),
            True,
            mult,
        )
        return RootIsolationResult(
            (exact_root,),
            complete=True,
            method="certified_interval_resolution",
            certificates=(cert,),
        )
    cert = RootCertificate(
        domain,
        1,
        "interval_resolution_rational_isolation",
        (*proof, ("rational_bracket", a, b), ("isolating_interval", isolating)),
        True,
        1,
    )
    return RootIsolationResult(
        isolated_intervals=(CertifiedRootInterval(isolating, cert),),
        complete=True,
        method="certified_interval_resolution",
        certificates=(cert,),
    )


def _certified_interval_isolation(expr, var, domain):
    """Certify all real roots on a finite interval from a complete monotone partition.

    Candidate discovery is delegated to the existing exact derivative-zero/monotonicity
    machinery.  Certification itself uses exact endpoint signs and rational
    opposite-sign brackets.  Completeness is asserted only when every monotone cell
    is classified and the partition is complete.
    """
    if (
        not isinstance(domain, sp.Interval)
        or domain.start in (-sp.oo, sp.oo)
        or domain.end in (-sp.oo, sp.oo)
    ):
        return None
    try:
        from ..real_analysis.monotonic_certificates import monotonic_partition

        partition = monotonic_partition(expr, var, domain=domain, return_result=True)
    except (ImportError, *SYMBOLIC_QUERY_ERRORS):
        return None
    if not partition.finite or not partition.cells:
        return None

    roots = []
    isolated = []
    certs = []
    complete = bool(partition.complete)

    # Critical points and attained outer endpoints can themselves be roots,
    # including even-multiplicity roots that do not produce a sign change.
    boundary_points = list(partition.critical_points)
    if not domain.left_open:
        boundary_points.append(domain.start)
    if not domain.right_open:
        boundary_points.append(domain.end)
    for point in _sort_unique(boundary_points):
        try:
            if sp.Contains(point, domain) is S.false:
                continue
            value = expr.subs(var, point)
        except SYMBOLIC_QUERY_ERRORS:
            complete = False
            continue
        if _exact_sign(value) == 0:
            roots.append(point)
            mult = _certified_multiplicity(expr, var, point)
            certs.append(
                RootCertificate(
                    sp.FiniteSet(point),
                    1,
                    "exact_partition_boundary_root",
                    (("value", value), ("partition", partition.source)),
                    True,
                    mult,
                )
            )

    for cell in partition.cells:
        interval = cell.interval
        if (
            not isinstance(interval, sp.Interval)
            or not cell.strict
            or not cell.continuous
        ):
            complete = False
            continue
        sl = _exact_sign(cell.left_limit) if cell.left_limit is not None else None
        sr = _exact_sign(cell.right_limit) if cell.right_limit is not None else None
        if sl is None or sr is None:
            complete = False
            continue
        # A zero limit at an excluded critical boundary is handled above if the
        # boundary belongs to the requested domain; strict monotonicity excludes
        # another zero in the open cell on the same side.
        if sl == 0 or sr == 0:
            certs.append(
                RootCertificate(
                    interval,
                    0,
                    "strict_monotone_boundary_limit",
                    (("endpoint_signs", sl, sr), ("direction", cell.direction)),
                    True,
                )
            )
            continue
        if sl == sr:
            certs.append(
                RootCertificate(
                    interval,
                    0,
                    "strict_monotone_same_sign",
                    (("endpoint_signs", sl, sr), ("direction", cell.direction)),
                    True,
                )
            )
            continue
        bracket = _rational_sign_bracket(expr, var, interval, sl, sr)
        if bracket is None:
            complete = False
            continue
        a, b = bracket
        exact_root, isolating = _bisect_certified_bracket(expr, var, a, b, sl)
        proof = (
            ("endpoint_signs", sl, sr),
            ("rational_bracket", a, b),
            ("direction", cell.direction),
            ("partition", partition.source),
        )
        if exact_root is not None:
            roots.append(exact_root)
            mult = _certified_multiplicity(expr, var, exact_root)
            certs.append(
                RootCertificate(
                    interval,
                    1,
                    "strict_monotonicity+exact_bisection",
                    (*proof, ("exact_root", exact_root)),
                    True,
                    mult,
                )
            )
        else:
            cert = RootCertificate(
                interval,
                1,
                "strict_monotonicity+rational_isolation",
                (*proof, ("isolating_interval", isolating)),
                True,
                1,
            )
            isolated.append(CertifiedRootInterval(isolating, cert))
            certs.append(cert)

    roots = _sort_unique(roots)
    method = (
        "certified_interval_isolation"
        if complete
        else "partial_certified_interval_isolation"
    )
    return RootIsolationResult(
        roots,
        isolated_intervals=tuple(isolated),
        complete=complete,
        method=method,
        evidence=(partition,),
        certificates=tuple(certs),
    )


def _scale_root_certificate(cert, factor, *, method_prefix="repeated_factor"):
    """Lift a certificate through ``g(x)**factor`` without losing provenance."""
    base_multiplicity = cert.multiplicity
    if (
        base_multiplicity is None
        and cert.root_count == 1
        and cert.method
        in {
            "strict_monotonicity+isolation",
            "strict_monotonicity+rational_isolation",
            "strict_monotonicity+exact_bisection",
        }
    ):
        base_multiplicity = 1
    multiplicity = base_multiplicity * factor if base_multiplicity is not None else None
    return RootCertificate(
        cert.domain,
        cert.root_count,
        f"{method_prefix}:{cert.method}",
        (*cert.proof, ("multiplicity_scale", factor)),
        cert.complete,
        multiplicity,
    )


def _repeated_power_roots(expr, var, domain):
    """Lift certified roots of a repeated factor ``base**m``.

    This is especially important for even powers: the original expression has no
    sign change at a simple zero of ``base``, while ``base`` itself does.  The
    recursive query disables heuristic numerical roots; every lifted root therefore
    retains a proof from the underlying factor.
    """
    if not isinstance(expr, sp.Pow):
        return None
    base, exponent = expr.as_base_exp()
    if exponent.is_Integer is not True or exponent < 2 or base == expr:
        return None
    underlying = transcendental_roots(base, var, domain, allow_numerical=False)
    if not (
        underlying.roots
        or underlying.isolated_intervals
        or underlying.zero_intervals
        or underlying.periodic is not None
        or underlying.complete
    ):
        return None
    factor = int(exponent)
    lifted_isolated = tuple(
        CertifiedRootInterval(
            item.interval, _scale_root_certificate(item.certificate, factor)
        )
        for item in underlying.isolated_intervals
    )
    lifted_certs = tuple(
        _scale_root_certificate(cert, factor) for cert in underlying.certificates
    )
    # Exact roots can often prove a sharper multiplicity directly (for example
    # when the base root itself is multiple), so prefer that over mechanical scaling.
    exact_certs = []
    for root in underlying.roots:
        multiplicity = _certified_multiplicity(expr, var, root)
        exact_certs.append(
            RootCertificate(
                sp.FiniteSet(root),
                1,
                "repeated_factor_exact_root",
                (("base", base), ("power", factor)),
                True,
                multiplicity,
            )
        )
    return RootIsolationResult(
        underlying.roots,
        underlying.zero_intervals,
        lifted_isolated,
        underlying.periodic,
        underlying.complete,
        "repeated_factor_root_isolation",
        (*underlying.evidence, ("repeated_factor", base, factor)),
        tuple(exact_certs) + lifted_certs,
        underlying.parametric_families,
        underlying.geometry,
    )


def _validated_even_multiplicity_numeric(expr, var, domain, *, max_samples=96):
    """Discover touching roots numerically and promote only exact validations.

    Numerical data is used solely to propose critical-point candidates.  A candidate
    becomes a certified root only when exact substitution proves ``f(r) == 0`` and
    repeated differentiation proves an even multiplicity.  Consequently this helper
    can improve coverage without treating a small residual as mathematical equality.
    """
    if (
        not isinstance(domain, sp.Interval)
        or domain.start in (-sp.oo, sp.oo)
        or domain.end in (-sp.oo, sp.oo)
    ):
        return None
    try:
        import mpmath as mp

        a = mp.mpf(str(sp.N(domain.start, 50)))
        b = mp.mpf(str(sp.N(domain.end, 50)))
        if not (mp.isfinite(a) and mp.isfinite(b) and a < b):
            return None
        derivative = sp.diff(expr, var)
        f = sp.lambdify(var, expr, "mpmath")
        df = sp.lambdify(var, derivative, "mpmath")
    except (ImportError, *SYMBOLIC_QUERY_ERRORS):
        return None

    old_dps = mp.mp.dps
    mp.mp.dps = max(old_dps, 80)
    try:
        xs = [a + (b - a) * i / max_samples for i in range(max_samples + 1)]
        candidates = []
        dvals = []
        for q in xs:
            try:
                value = df(q)
                dvals.append(mp.mpf(value) if mp.im(value) == 0 else mp.nan)
            except (ArithmeticError, ValueError, TypeError):
                dvals.append(mp.nan)
        for x1, x2, y1, y2 in zip(xs, xs[1:], dvals, dvals[1:], strict=False):
            if not (mp.isfinite(y1) and mp.isfinite(y2)):
                continue
            if y1 == 0:
                candidates.append(x1)
            if y1 * y2 < 0:
                try:
                    candidates.append(mp.findroot(df, (x1, x2), solver="anderson"))
                except (ArithmeticError, ValueError, ZeroDivisionError):
                    with contextlib.suppress(
                        ArithmeticError, ValueError, ZeroDivisionError
                    ):
                        candidates.append(mp.findroot(df, (x1, x2)))
        # Sample-local minima of |f| catch stationary roots whose derivative root is
        # numerically flat enough that a sign-change scan of f' misses it.
        fvals = []
        for q in xs:
            try:
                value = f(q)
                fvals.append(abs(value) if mp.im(value) == 0 else mp.inf)
            except (ArithmeticError, ValueError, TypeError):
                fvals.append(mp.inf)
        scale = max([v for v in fvals if mp.isfinite(v)] or [mp.mpf(1)])
        threshold = max(mp.mpf("1e-30"), scale * mp.mpf("1e-18"))
        for i in range(1, len(xs) - 1):
            if (
                fvals[i] <= fvals[i - 1]
                and fvals[i] <= fvals[i + 1]
                and fvals[i] <= threshold
            ):
                with contextlib.suppress(
                    ArithmeticError, ValueError, ZeroDivisionError
                ):
                    candidates.append(mp.findroot(df, (xs[i - 1], xs[i + 1])))

        roots = []
        certs = []
        for candidate in candidates:
            if not mp.isfinite(candidate) or not (a <= candidate <= b):
                continue
            approx = sp.Float(str(candidate), 70)
            guesses = []
            for constants in ((), (sp.pi, sp.E, sp.sqrt(2), sp.sqrt(3), sp.log(2))):
                with contextlib.suppress(*SYMBOLIC_QUERY_ERRORS):
                    guesses.append(
                        sp.nsimplify(
                            approx,
                            constants=list(constants),
                            tolerance=1e-28,
                            full=bool(constants),
                        )
                    )
            # Keep the exact decimal rational as a final candidate.  It almost never
            # validates, but makes the reconstruction policy explicit and deterministic.
            guesses.append(sp.Rational(str(candidate)))
            for root in guesses:
                try:
                    if sp.Contains(root, domain) is S.false:
                        continue
                    if _exact_sign(expr.subs(var, root)) != 0:
                        continue
                    multiplicity = _certified_multiplicity(
                        expr, var, root, max_order=16
                    )
                except SYMBOLIC_QUERY_ERRORS:
                    continue
                if multiplicity is None or multiplicity < 2 or multiplicity % 2:
                    continue
                if root not in roots:
                    roots.append(root)
                    certs.append(
                        RootCertificate(
                            sp.FiniteSet(root),
                            1,
                            "numerical_candidate+exact_even_multiplicity_validation",
                            (
                                ("approximation", approx),
                                ("derivative_candidate", derivative),
                            ),
                            True,
                            multiplicity,
                        )
                    )
                break
        if not roots:
            return None
        return RootIsolationResult(
            _sort_unique(roots),
            complete=False,
            method="validated_even_multiplicity_numeric_fallback",
            evidence=(("numeric_candidates", len(candidates)),),
            certificates=tuple(certs),
        )
    finally:
        mp.mp.dps = old_dps


def _numeric(expr, var, domain, max_samples=80):
    if (
        not isinstance(domain, sp.Interval)
        or domain.start in (-sp.oo, sp.oo)
        or domain.end in (-sp.oo, sp.oo)
    ):
        return None
    try:
        a = float(sp.N(domain.start))
        b = float(sp.N(domain.end))
    except SYMBOLIC_QUERY_ERRORS:
        return None
    if not math.isfinite(a) or not math.isfinite(b) or not a < b:
        return None
    import mpmath as mp

    f = sp.lambdify(var, expr, "mpmath")
    roots = []
    xs = [a + (b - a) * i / max_samples for i in range(max_samples + 1)]
    vals = []
    for q in xs:
        try:
            vals.append(float(f(q)))
        except SYMBOLIC_QUERY_ERRORS:
            vals.append(float("nan"))
    for x1, x2, y1, y2 in zip(xs, xs[1:], vals, vals[1:], strict=False):
        if not (math.isfinite(y1) and math.isfinite(y2)):
            continue
        if y1 * y2 < 0:
            with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
                roots.append(sp.Float(str(mp.findroot(lambda t: f(t), (x1, x2))), 30))
    return RootIsolationResult(
        _sort_unique(roots),
        complete=False,
        method="numerical_scan",
        evidence=("candidate roots only",),
    )


def _exponential_fourier_count(expr, var, domain):
    """Attach a certified derivative-sequence count for supported exp polynomials."""
    try:
        from .exp_fourier import exponential_fourier_count

        count = exponential_fourier_count(expr, var, domain)
    except (ImportError, *SYMBOLIC_QUERY_ERRORS):
        return None
    if count is None or not count.complete or count.root_count is None:
        return None
    cert = RootCertificate(
        domain,
        count.root_count,
        "exponential_fourier_derivative_sequence",
        (count,),
        True,
    )
    if count.root_count == 0:
        return RootIsolationResult(
            complete=True,
            method="exponential_fourier_count",
            evidence=(count,),
            certificates=(cert,),
        )
    return RootIsolationResult(
        complete=False,
        method="exponential_fourier_count",
        evidence=(count,),
        certificates=(cert,),
    )


def transcendental_roots(expr, var, domain=S.Reals, *, allow_numerical=True):
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if expr == 0:
        zi = (domain,) if isinstance(domain, sp.Interval) else ()
        c = RootCertificate(domain, None, "identically_zero", (expr,), True)
        return RootIsolationResult(
            (), zi, complete=True, method="identically_zero", certificates=(c,)
        )
    if not expr.has(var):
        c = RootCertificate(domain, 0, "constant", (expr,), True)
        return RootIsolationResult(
            (), complete=True, method="constant", certificates=(c,)
        )
    # General zero-geometry theorems precede generic solveset.  They preserve
    # indexed families/interlacing even when individual zeros have no closed form.
    geom = zero_geometry(expr, var, domain=domain)
    if geom:
        exact_guard = all(_guard_proved(g.guard) for g in geom)
        support_ok = all(
            domain == g.support
            or (domain is S.Reals and g.family in ("AiryAi", "AiryBi"))
            for g in geom
        )
        priority = all(
            g.family in ("AiryAi", "AiryBi", "BesselJ", "BesselY", "FerrersQ")
            for g in geom
        )
        if priority and exact_guard and support_ok:
            certs = tuple(
                RootCertificate(
                    domain,
                    None,
                    "zero_geometry",
                    ((g.source, g.count, g.support),),
                    g.complete_on_support,
                )
                for g in geom
            )
            return RootIsolationResult(
                complete=all(g.complete_on_support for g in geom),
                method="zero_geometry",
                evidence=geom,
                certificates=certs,
                geometry=geom,
            )
    for fn in (
        _periodic,
        _special,
        _orthogonal_special,
        _repeated_power_roots,
        _identically_zero_interval,
        _solveset,
        _exponential_fourier_count,
        _certified_interval_resolution,
        _certified_monotone_interval,
        _certified_interval_isolation,
    ):
        r = fn(expr, var, domain)
        if r is not None:
            if geom and not r.geometry:
                r = RootIsolationResult(
                    r.roots,
                    r.zero_intervals,
                    r.isolated_intervals,
                    r.periodic,
                    r.complete,
                    r.method,
                    r.evidence,
                    r.certificates,
                    r.parametric_families,
                    geom,
                )
            return r
    if geom:
        return RootIsolationResult(
            complete=False,
            method="conditional_zero_geometry",
            evidence=geom,
            geometry=geom,
        )
    if allow_numerical:
        r = _validated_even_multiplicity_numeric(expr, var, domain)
        if r is not None:
            return r
        r = _numeric(expr, var, domain)
        if r is not None:
            return r
    return RootIsolationResult(complete=False, method="unknown")
