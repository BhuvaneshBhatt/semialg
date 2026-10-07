"""Certified definite and complex-contour integration.

The public :func:`integrate` returns the mathematical value by default.
``return_result=True`` exposes the replayable proof/result object.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import mpmath as mp
import sympy as sp
from sympy.printing.codeprinter import PrintMethodNotImplementedError

from .exactness import certify_identity_equal


@dataclass(frozen=True)
class ParametricContour:
    """Oriented contour ``gamma(t)``, ``a <= t <= b``.

    The parameterization is retained: reversing the parameter interval reverses
    orientation even when the geometric image is unchanged.
    """

    gamma: sp.Expr
    parameter: sp.Symbol
    lower: sp.Expr
    upper: sp.Expr

    def __init__(self, gamma, parameter_range):
        if len(parameter_range) != 3:
            raise ValueError("parameter_range must be (parameter, lower, upper)")
        t, a, b = parameter_range
        object.__setattr__(self, "gamma", sp.sympify(gamma))
        object.__setattr__(self, "parameter", sp.sympify(t))
        object.__setattr__(self, "lower", sp.sympify(a))
        object.__setattr__(self, "upper", sp.sympify(b))
        if not isinstance(self.parameter, sp.Symbol):
            raise TypeError("contour parameter must be a Symbol")
        if (
            self.parameter not in self.gamma.free_symbols
            and self.gamma.is_number is not True
        ):
            raise ValueError("gamma must depend on the contour parameter")

    @property
    def start(self):
        return sp.simplify(self.gamma.subs(self.parameter, self.lower))

    @property
    def end(self):
        return sp.simplify(self.gamma.subs(self.parameter, self.upper))

    @property
    def closed(self) -> bool:
        return certify_identity_equal(self.start, self.end) is True

    @property
    def derivative(self):
        return sp.diff(self.gamma, self.parameter)

    def reversed(self):
        return ParametricContour(self.gamma, (self.parameter, self.upper, self.lower))


@dataclass(frozen=True)
class IntegralResult:
    expression: sp.Expr
    variable: sp.Symbol
    path: object
    value: sp.Expr | None
    method: str
    complete: bool
    certified: bool
    error_bound: sp.Expr | None = None
    evidence: tuple = ()

    def verify(self) -> bool:
        if (
            self.method == "exact-definite"
            and isinstance(self.path, tuple)
            and len(self.path) == 2
        ):
            lim = (self.variable, self.path[0], self.path[1])
        else:
            lim = (self.variable, self.path)
        fresh = integrate(self.expression, lim, return_result=True, _replay=True)
        if not isinstance(fresh, IntegralResult):
            return False
        if fresh.complete != self.complete or fresh.certified != self.certified:
            return False
        if fresh.method != self.method:
            return False
        if self.value is None or fresh.value is None:
            return self.value is fresh.value
        return certify_identity_equal(fresh.value, self.value) is True or (
            self.error_bound is not None
            and fresh.error_bound == self.error_bound
            and abs(complex(sp.N(fresh.value - self.value, 30)))
            <= float(sp.N(self.error_bound, 20))
        )


def _unknown(expr, var, path, method, evidence=()):
    return IntegralResult(
        expr, var, path, None, method, False, False, None, tuple(evidence)
    )


def _exact_parameter_integral(expr, var, contour):
    t = contour.parameter
    pulled = sp.simplify(expr.subs(var, contour.gamma) * contour.derivative)
    # Risch is a fast, terminating first attempt for elementary pulled-backs.
    # It avoids letting a difficult exact definite integral monopolize the
    # certified numerical fallback.
    try:
        anti = sp.integrate(pulled, t, risch=True)
    except (ValueError, NotImplementedError):
        anti = sp.Integral(pulled, t)
    if not isinstance(anti, sp.Integral) and not anti.has(sp.Integral):
        # Endpoint subtraction of a principal antiderivative is unsound on a
        # closed complex path when that antiderivative is multivalued: it can
        # erase logarithmic/power monodromy. Keep such cases for a contour-
        # aware certificate.
        multivalued = anti.has(sp.log) or any(
            p.exp.is_integer is not True for p in anti.atoms(sp.Pow)
        )
        if not (contour.closed and multivalued):
            val = sp.simplify(anti.subs(t, contour.upper) - anti.subs(t, contour.lower))
            return val, pulled
    # Rational and algebraizable-trig cases have reliable exact definite paths.
    if _rational_parameterization(contour.gamma, t) or _trig_algebraizable(
        contour.gamma, t
    ):
        val = sp.integrate(pulled, (t, contour.lower, contour.upper))
        if not isinstance(val, sp.Integral) and not val.has(sp.Integral):
            return sp.simplify(val), pulled
    return None, pulled


def _rational_parameterization(gamma, t):
    re, im = sp.expand_complex(gamma).as_real_imag()
    try:
        return (
            sp.Poly(sp.together(re).as_numer_denom()[0], t),
            sp.Poly(sp.together(re).as_numer_denom()[1], t),
            sp.Poly(sp.together(im).as_numer_denom()[0], t),
            sp.Poly(sp.together(im).as_numer_denom()[1], t),
        )
    except sp.PolynomialError:
        return None


def _trig_algebraizable(gamma, t):
    allowed = {sp.sin, sp.cos}
    for f in gamma.atoms(sp.Function):
        if f.func not in allowed:
            return False
        arg = sp.expand(f.args[0])
        if sp.diff(arg, t, 2) != 0:
            return False
    return bool(gamma.atoms(sp.sin, sp.cos))


def _trig_algebraization(gamma, t):
    """Return an algebraic unit-circle representation when trig arguments are integer-affine."""
    if not _trig_algebraizable(gamma, t):
        return None
    X, Y = sp.symbols("_contour_cos _contour_sin", real=True)
    expanded = sp.expand_trig(gamma)
    # expand_trig reduces integer multiples/shifts to sin(t), cos(t) plus constants.
    alg = sp.expand(expanded.xreplace({sp.sin(t): Y, sp.cos(t): X}))
    if alg.has(sp.sin(t), sp.cos(t)):
        return None
    return alg, sp.Eq(X**2 + Y**2, 1), X, Y


def _point_avoided_semialg(contour, point):
    """Prove a rational contour avoids a point; delegate satisfiability to semialg."""
    t = contour.parameter
    if _rational_parameterization(contour.gamma, t) is None:
        return None
    re, im = sp.expand_complex(contour.gamma - point).as_real_imag()
    nr, dr = map(sp.expand, sp.together(re).as_numer_denom())
    ni, di = map(sp.expand, sp.together(im).as_numer_denom())
    interval = sp.And(
        t >= sp.Min(contour.lower, contour.upper),
        t <= sp.Max(contour.lower, contour.upper),
    )
    cond = sp.And(sp.Eq(nr, 0), sp.Eq(ni, 0), sp.Ne(dr, 0), sp.Ne(di, 0), interval)
    try:
        from ._semialg import is_satisfiable

        return not is_satisfiable(cond, (t,))
    except (ImportError, AttributeError, TypeError, ValueError, NotImplementedError):
        return None


def _winding_exact(contour, point):
    if not contour.closed:
        return None
    t = contour.parameter
    q = sp.simplify(contour.derivative / (contour.gamma - point))
    v = sp.integrate(q, (t, contour.lower, contour.upper))
    if isinstance(v, sp.Integral) or v.has(sp.Integral):
        return None
    n = sp.simplify(v / (2 * sp.pi * sp.I))
    if n.is_integer is True:
        return int(n)
    if n.is_Integer:
        return int(n)
    return None


def _interval_derivative_bound(expr, t, lo, hi):
    """Validated mpmath interval upper bound for |expr| on a real panel."""
    try:
        f = sp.lambdify(t, expr, modules=mp.iv)
        lo_f = math.nextafter(float(sp.N(lo, 30)), -math.inf)
        hi_f = math.nextafter(float(sp.N(hi, 30)), math.inf)
        v = abs(f(mp.iv.mpf([lo_f, hi_f])))
        upper = float(v.b)
        # One binary64 ulp outward preserves the interval upper endpoint after conversion.
        return sp.Float(math.nextafter(upper, math.inf), 17)
    except (
        ValueError,
        TypeError,
        AttributeError,
        ZeroDivisionError,
        OverflowError,
        PrintMethodNotImplementedError,
    ):
        return None


def _adaptive_taylor_integral(pulled, t, a, b, *, tol=None, max_depth=10, order=6):
    """Certified Taylor-panel quadrature when SymPy can bound the remainder.

    Each panel expands at its midpoint.  A global maximum of |h^(n+1)| on the
    real panel supplies the Lagrange remainder bound. Unsupported maxima cause
    conservative failure.
    """
    if tol is None:
        tol = sp.Rational(1, 10) ** 12
    tol = sp.sympify(tol)
    derivatives = [sp.diff(pulled, t, k) for k in range(order + 2)]

    def panel(lo, hi, depth):
        mid = sp.simplify((lo + hi) / 2)
        rad = sp.simplify(abs(hi - lo) / 2)
        val = sp.S.Zero
        for k in range(order + 1):
            ck = derivatives[k].subs(t, mid) / sp.factorial(k)
            if ck.has(sp.zoo, sp.nan, sp.oo, -sp.oo):
                return None
            val += (
                ck * ((hi - mid) ** (k + 1) - (lo - mid) ** (k + 1)) / sp.Integer(k + 1)
            )
        d = derivatives[order + 1]
        # Prefer validated interval arithmetic; symbolic maxima are a fallback.
        M = _interval_derivative_bound(d, t, sp.Min(lo, hi), sp.Max(lo, hi))
        if M is None:
            return None
        err = sp.simplify(2 * M * rad ** (order + 2) / sp.factorial(order + 2))
        if err.is_finite is not True or err.is_nonnegative is not True:
            # SymPy may fail to decide Float finiteness; numeric finite nonnegative is enough here.
            try:
                if not (float(sp.N(err)) >= 0 and float(sp.N(err)) < float("inf")):
                    return None
            except (TypeError, ValueError, OverflowError):
                return None
        try:
            accepted = float(sp.N(err, 20)) <= float(sp.N(tol, 20))
        except (TypeError, ValueError, OverflowError):
            accepted = False
        if accepted:
            return val, err, 1
        if depth >= max_depth:
            return None
        left = panel(lo, mid, depth + 1)
        right = panel(mid, hi, depth + 1)
        if left is None or right is None:
            return None
        return left[0] + right[0], left[1] + right[1], left[2] + right[2]

    return panel(sp.sympify(a), sp.sympify(b), 0)


def _rational_path_clear(expr, var, contour):
    """Prove the pulled-back rational integrand has no denominator zero."""
    t = contour.parameter
    pulled = sp.cancel(expr.subs(var, contour.gamma) * contour.derivative)
    num, den = sp.fraction(pulled)
    try:
        sp.Poly(num, t)
        sp.Poly(den, t)
    except sp.PolynomialError:
        return None
    if sp.degree(den, t) <= 0:
        return True
    lo, hi = sp.Min(contour.lower, contour.upper), sp.Max(contour.lower, contour.upper)
    cond = sp.And(sp.Eq(den, 0), t >= lo, t <= hi)
    try:
        from ._semialg import is_satisfiable

        return not is_satisfiable(cond, (t,))
    except (ImportError, AttributeError, TypeError, ValueError, NotImplementedError):
        roots = sp.solveset(sp.Eq(den, 0), t, domain=sp.Interval(lo, hi))
        return roots is sp.S.EmptySet


def _branch_sensitive(expr):
    if expr.has(sp.log):
        return True
    for p in expr.atoms(sp.Pow):
        if p.exp.is_integer is not True:
            return True
    return False


def _simple_positive_real_path(contour):
    t = contour.parameter
    im = sp.simplify(sp.im(sp.expand_complex(contour.gamma)))
    if im != 0:
        return False
    # Sufficient endpoint/monotone-affine case; otherwise remain conservative.
    try:
        q = sp.Poly(sp.re(sp.expand_complex(contour.gamma)), t)
        if q.degree() > 1:
            return False
        vals = [
            sp.simplify(q.as_expr().subs(t, u)) for u in (contour.lower, contour.upper)
        ]
        return all(v.is_positive is True for v in vals)
    except sp.PolynomialError:
        return False


def _pulled_continuous_on_interval(pulled, t, a, b):
    """Prove the pulled-back integrand is continuous on the closed panel."""
    interval = sp.Interval(sp.Min(a, b), sp.Max(a, b))
    try:
        from sympy.calculus.util import continuous_domain

        domain = continuous_domain(pulled, t, interval)
    except (TypeError, ValueError, NotImplementedError, ImportError, RuntimeError):
        return None
    return domain == interval


def _parametric(expr, var, contour, *, tolerance=None):
    exact, pulled = _exact_parameter_integral(expr, var, contour)
    if _branch_sensitive(expr) and not _simple_positive_real_path(contour):
        return _unknown(
            expr,
            var,
            contour,
            "branch-sensitive-path-unresolved",
            (("pulled_back", pulled),),
        )
    if exact is not None:
        if exact.has(sp.nan, sp.zoo, sp.oo, -sp.oo) or exact.is_finite is False:
            return _unknown(
                expr,
                var,
                contour,
                "parametric-singular-path",
                (("pulled_back", pulled),),
            )
        if (
            _pulled_continuous_on_interval(
                pulled, contour.parameter, contour.lower, contour.upper
            )
            is not True
        ):
            return _unknown(
                expr,
                var,
                contour,
                "parametric-path-continuity-unproved",
                (("pulled_back", pulled),),
            )
        rational = _rational_parameterization(contour.gamma, contour.parameter)
        if (
            rational
            and sp.cancel(pulled).is_rational_function(contour.parameter)
            and _rational_path_clear(expr, var, contour) is not True
        ):
            return _unknown(
                expr,
                var,
                contour,
                "parametric-rational-path-unclear",
                (("pulled_back", pulled),),
            )
        trig = (
            _trig_algebraization(contour.gamma, contour.parameter)
            if not rational
            else None
        )
        kind = (
            "parametric-rational-exact"
            if rational
            else (
                "parametric-trig-algebraized-exact"
                if trig is not None
                else "parametric-analytic-exact"
            )
        )
        evidence = [("pulled_back", pulled)]
        if trig is not None:
            evidence.append(("trig_algebraization", trig[:2]))
        return IntegralResult(
            expr, var, contour, exact, kind, True, True, evidence=tuple(evidence)
        )
    tol = sp.Rational(1, 10) ** 12 if tolerance is None else sp.sympify(tolerance)
    q = _adaptive_taylor_integral(
        pulled, contour.parameter, contour.lower, contour.upper, tol=tol
    )
    if q is None:
        return _unknown(
            expr, var, contour, "parametric-analytic-taylor", (("pulled_back", pulled),)
        )
    val, err, panels = q
    return IntegralResult(
        expr,
        var,
        contour,
        sp.N(val, 30),
        "parametric-analytic-taylor",
        True,
        True,
        err,
        (("pulled_back", pulled), ("panels", panels)),
    )


def _polygon_contour(points):
    pts = tuple(map(sp.sympify, points))
    if len(pts) < 2:
        raise ValueError("a polygonal path needs at least two points")
    return pts


def _polygon_integral(expr, var, points):
    pts = _polygon_contour(points)
    total = sp.S.Zero
    for i, (a, b) in enumerate(itertools.pairwise(pts)):
        t = sp.Dummy(f"_s{i}", real=True)
        gamma = a + t * (b - a)
        c = ParametricContour(gamma, (t, 0, 1))
        rr = _parametric(expr, var, c)
        if not rr.complete or rr.value is None:
            return _unknown(expr, var, pts, "polygon-segment-exact")
        total = sp.simplify(total + rr.value)
    return IntegralResult(expr, var, pts, total, "polygon-segment-exact", True, True)


def integrate(expr, limits, *, return_result=False, tolerance=None, _replay=False):
    """Integrate with conservative certification.

    ``integrate(expr, (x,a,b))`` integrates over a finite interval/straight
    complex segment. ``integrate(expr, (z, contour))`` accepts a
    :class:`ParametricContour` or a sequence of polygon vertices.
    """
    expr = sp.sympify(expr)
    if not isinstance(limits, (tuple, list)) or len(limits) not in (2, 3):
        raise TypeError("limits must be (variable, path) or (variable, lower, upper)")
    var = sp.sympify(limits[0])
    if not isinstance(var, sp.Symbol):
        raise TypeError("integration variable must be a Symbol")
    if len(limits) == 3:
        a, b = map(sp.sympify, limits[1:])
        t = sp.Dummy("_segment_t", real=True)
        contour = ParametricContour(a + t * (b - a), (t, 0, 1))
        result = _parametric(expr, var, contour, tolerance=tolerance)
        # Preserve conventional exact method name for real/complex straight limits.
        if result.complete and result.method.endswith("exact"):
            result = IntegralResult(
                expr,
                var,
                (a, b),
                result.value,
                "exact-definite",
                True,
                True,
                evidence=result.evidence,
            )
    else:
        path = limits[1]
        if isinstance(path, ParametricContour):
            result = _parametric(expr, var, path, tolerance=tolerance)
        elif isinstance(path, (tuple, list)):
            result = _polygon_integral(expr, var, path)
        else:
            raise TypeError(
                "path must be ParametricContour or a sequence of complex points"
            )
    if return_result:
        return result
    return (
        result.value
        if result.complete
        else sp.Integral(expr, limits if len(limits) == 3 else (var,))
    )
