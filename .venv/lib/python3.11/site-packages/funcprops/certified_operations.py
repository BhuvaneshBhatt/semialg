"""Replayable high-level certified operations built from Funcprops facts."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .exactness import certify_identity_equal
from .integration import ParametricContour, integrate


@dataclass(frozen=True)
class InverseCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    value: sp.Symbol
    inverse: sp.Expr
    domain: sp.Set
    image: sp.Set

    def verify(self):
        from sympy.calculus.util import function_range as sympy_function_range

        if not isinstance(self.value, sp.Symbol):
            return False
        try:
            actual_image = sympy_function_range(
                self.expression, self.variable, self.domain
            )
        except (TypeError, ValueError, NotImplementedError, ImportError, RuntimeError):
            return False
        if actual_image != self.image:
            return False
        right = (
            certify_identity_equal(
                self.expression.subs(self.variable, self.inverse), self.value
            )
            is True
        )
        left = (
            certify_identity_equal(
                self.inverse.subs(self.value, self.expression), self.variable
            )
            is True
        )
        return right and left


def inverse(expr, var, *, value=None, domain=sp.S.Reals, return_certificate=False):
    e = sp.sympify(expr)
    x = sp.sympify(var)
    y = sp.Dummy("_inverse_value") if value is None else sp.sympify(value)
    if not isinstance(y, sp.Symbol):
        raise TypeError("value must be a Symbol naming the inverse-function argument")
    from .api import injective

    inj = injective(e, x, domain=domain)
    # A unique symbolic inverse plus a certified nonzero constant derivative is
    # itself an injectivity proof, avoiding dependence on result presentation.
    d = sp.simplify(sp.diff(e, x))
    affine_injective = not d.has(x) and d.is_zero is False
    sol = sp.solve(sp.Eq(e, y), x)
    inj_true = (
        inj is True
        or getattr(inj, "name", "") in ("TRUE", "INJECTIVE")
        or affine_injective
    )
    if len(sol) != 1 or not inj_true:
        return None
    inv = sp.simplify(sol[0])
    # The certificate stores the actual image, so do not substitute a guessed
    # UniversalSet when range computation is unavailable.
    try:
        from sympy.calculus.util import function_range as sympy_function_range

        img = sympy_function_range(e, x, domain)
    except (TypeError, ValueError, NotImplementedError, ImportError, RuntimeError):
        return None
    if not isinstance(img, sp.Set):
        return None
    c = InverseCertificate(e, x, y, inv, domain, img)
    if not c.verify():
        return None
    return c if return_certificate else inv


@dataclass(frozen=True)
class EquationSolutionCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    domain: sp.Set
    solutions: sp.Set
    complete: bool

    def verify(self):
        s = sp.solveset(self.expression, self.variable, domain=self.domain)
        complete = not isinstance(s, sp.ConditionSet)
        return s == self.solutions and self.complete == complete


def solve(expr, var, *, domain=sp.S.Complexes, return_certificate=False):
    e = sp.sympify(expr)
    x = sp.sympify(var)
    s = sp.solveset(e, x, domain=domain)
    complete = not isinstance(s, sp.ConditionSet)
    c = EquationSolutionCertificate(e, x, domain, s, complete)
    return c if return_certificate else s


@dataclass(frozen=True)
class LaurentCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    point: sp.Expr
    series: sp.Expr
    order: int

    def verify(self):
        return (
            sp.series(self.expression, self.variable, self.point, self.order)
            == self.series
        )


def laurent_series(expr, var, point=0, *, order=6, return_certificate=False):
    e = sp.sympify(expr)
    x = sp.sympify(var)
    p = sp.sympify(point)
    s = sp.series(e, x, p, order)
    c = LaurentCertificate(e, x, p, s, order)
    return c if return_certificate else s


@dataclass(frozen=True)
class ResidueCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    point: sp.Expr
    residue: sp.Expr

    def verify(self):
        return (
            certify_identity_equal(
                sp.residue(self.expression, self.variable, self.point), self.residue
            )
            is True
        )


def residue(expr, var, point, *, return_certificate=False):
    e = sp.sympify(expr)
    x = sp.sympify(var)
    p = sp.sympify(point)
    r = sp.residue(e, x, p)
    c = ResidueCertificate(e, x, p, r)
    return c if return_certificate else r


@dataclass(frozen=True)
class WindingCertificate:
    contour: ParametricContour
    point: sp.Expr
    winding: int
    integral: sp.Expr

    def verify(self):
        expected = 2 * sp.pi * sp.I * self.winding
        return (
            winding_number(self.contour, self.point) == self.winding
            and certify_identity_equal(self.integral, expected) is True
        )


def winding_number(contour, point=0, *, return_certificate=False):
    if not isinstance(contour, ParametricContour) or not contour.closed:
        return None
    p = sp.sympify(point)
    t = contour.parameter
    g = contour.gamma
    # Exact affine circle recognition c + A*exp(I*(n*t+phi)).  This path must
    # precede symbolic integration: definite log antiderivatives can erase turns.
    exps = [a for a in g.atoms(sp.exp) if a.args[0].has(t)]
    if len(exps) == 1:
        E = exps[0]
        arg = sp.expand(E.args[0])
        phi = sp.simplify(arg / sp.I)
        # exp(i*phi(t)) has winding (phi(b)-phi(a))/(2*pi) when phi is real
        # on the real parameter interval; no monotonicity is required.
        turns = sp.simplify(
            (phi.subs(t, contour.upper) - phi.subs(t, contour.lower)) / (2 * sp.pi)
        )
        if turns.is_Integer is True and (phi.is_real is True or not phi.has(sp.I)):
            n = int(turns)
            A = sp.simplify(g.coeff(E))
            center = sp.simplify(g - A * E)
            if not center.has(t) and not A.has(t):
                delta = sp.simplify(abs(A) ** 2 - sp.Abs(p - center) ** 2)
                if delta.is_positive is True:
                    w = n
                elif delta.is_negative is True:
                    w = 0
                else:
                    w = None
                if w is not None:
                    c = WindingCertificate(contour, p, w, 2 * sp.pi * sp.I * w)
                    return c if return_certificate else w
    z = sp.Dummy("_z")
    r = integrate(1 / (z - p), (z, contour), return_result=True)
    if not r.complete:
        return None
    n = sp.simplify(r.value / (2 * sp.pi * sp.I))
    if n.is_Integer is not True:
        return None
    c = WindingCertificate(contour, p, int(n), r.value)
    return c if return_certificate else c.winding


@dataclass(frozen=True)
class ArgumentPrincipleResult:
    expression: sp.Expr
    variable: sp.Symbol
    contour: ParametricContour
    zeros_minus_poles: int
    certified: bool

    def verify(self):
        fresh = argument_principle(self.expression, self.variable, self.contour)
        return self.certified is True and fresh == self.zeros_minus_poles


def argument_principle(expr, var, contour, *, return_certificate=False):
    e = sp.cancel(sp.sympify(expr))
    x = sp.sympify(var)
    # For rational meromorphic functions, sum oriented winding numbers of the
    # exact divisor. This avoids antiderivative branch-loss on closed contours.
    num, den = sp.fraction(e)
    total = 0
    try:
        pn, pd = sp.Poly(num, x), sp.Poly(den, x)
        nr, dr = sp.roots(pn), sp.roots(pd)
        if sum(nr.values()) != pn.degree() or sum(dr.values()) != pd.degree():
            return None
        for root, m in nr.items():
            w = winding_number(contour, root)
            if w is None:
                return None
            total += m * w
        for root, m in dr.items():
            w = winding_number(contour, root)
            if w is None:
                return None
            total -= m * w
        c = ArgumentPrincipleResult(e, x, contour, int(total), True)
        return c if return_certificate else c.zeros_minus_poles
    except sp.PolynomialError:
        # An integer value of integral(f'/f) is not by itself a proof of the
        # argument principle: meromorphicity in the enclosed region and
        # nonvanishing on the contour are theorem hypotheses.  The exact
        # rational divisor path above certifies those obligations structurally;
        # unsupported transcendental cases remain unknown.
        return None


@dataclass(frozen=True)
class HomotopyCertificate:
    """Certificate for based-loop homotopy in a finitely punctured plane.

    ``word`` is the signed crossing word for a certified system of pairwise
    disjoint parallel cuts from the punctures to infinity.  ``reduced_word``
    is its free reduction.  Equality of reduced words is the complete
    invariant for based loops in ``C \\ punctures``.
    """

    first: ParametricContour
    second: ParametricContour
    punctures: tuple
    cut_direction: sp.Expr | None
    first_word: tuple
    second_word: tuple
    reduced_word: tuple

    def verify(self):
        if (
            certify_identity_equal(self.first.start, self.second.start) is not True
            or certify_identity_equal(self.first.end, self.second.end) is not True
            or not self.first.closed
            or not self.second.closed
        ):
            return False
        ps = _distinct_punctures(self.punctures)
        if ps != self.punctures:
            return False
        if self.cut_direction is None:
            if len(ps) > 1:
                return False
            if not ps:
                a = b = ()
            else:
                wa = winding_number(self.first, ps[0])
                wb = winding_number(self.second, ps[0])
                if wa is None or wb is None:
                    return False
                a = tuple(([1] * wa) if wa >= 0 else ([-1] * (-wa)))
                b = tuple(([1] * wb) if wb >= 0 else ([-1] * (-wb)))
        else:
            if not _cut_system_valid(ps, self.first.start, self.cut_direction):
                return False
            a = _cut_crossing_word(self.first, ps, self.cut_direction)
            b = _cut_crossing_word(self.second, ps, self.cut_direction)
            if a is None or b is None:
                return False
        ar = _free_reduce(a)
        br = _free_reduce(b)
        return (
            a == self.first_word
            and b == self.second_word
            and ar == br == self.reduced_word
        )


def _distinct_punctures(punctures):
    out = []
    for p in map(sp.sympify, punctures):
        if not any(certify_identity_equal(p, q) is True for q in out):
            out.append(p)
    return tuple(out)


def _real_imag(expr):
    """Return exact real/imaginary parts after complex expansion."""
    return tuple(map(sp.simplify, sp.expand_complex(expr).as_real_imag()))


def _strict_sign(expr):
    e = sp.simplify(expr)
    if e.is_positive is True:
        return 1
    if e.is_negative is True:
        return -1
    if e.is_zero is True:
        return 0
    return None


def _cut_system_valid(punctures, basepoint, direction):
    """Prove that parallel cuts p + s*d (s>=0) are disjoint and miss basepoint."""
    d = sp.sympify(direction)
    if certify_identity_equal(d, 0) is True:
        return False
    dc = sp.conjugate(d)
    # Parallel rays are disjoint when no pair of punctures lies on a line
    # parallel to the cut direction.
    for i, p in enumerate(punctures):
        for q in punctures[i + 1 :]:
            _, transverse = _real_imag((p - q) * dc)
            if _strict_sign(transverse) in (None, 0):
                return False
        along, transverse = _real_imag((sp.sympify(basepoint) - p) * dc)
        ts = _strict_sign(transverse)
        if ts == 0 and _strict_sign(along) in (0, 1):
            return False
        if ts is None:
            return False
    return True


def _interior(root, lower, upper):
    """Prove root lies strictly inside the oriented parameter interval."""
    lo = sp.simplify(root - lower)
    hi = sp.simplify(upper - root)
    if _strict_sign(sp.simplify(upper - lower)) == 1:
        return _strict_sign(lo) == 1 and _strict_sign(hi) == 1
    if _strict_sign(sp.simplify(lower - upper)) == 1:
        return _strict_sign(-lo) == 1 and _strict_sign(-hi) == 1
    return False


def _ordered_events(events, lower, upper):
    """Exactly order (parameter, generator, sign) events along orientation."""
    reverse = _strict_sign(sp.simplify(upper - lower)) == -1
    pending = list(events)
    ordered = []
    while pending:
        chosen = None
        for i, event in enumerate(pending):
            t0 = event[0]
            ok = True
            for j, other in enumerate(pending):
                if i == j:
                    continue
                diff = sp.simplify(t0 - other[0])
                sg = _strict_sign(diff)
                if sg is None or sg == 0:
                    ok = False
                    break
                if (not reverse and sg > 0) or (reverse and sg < 0):
                    ok = False
                    break
            if ok:
                chosen = i
                break
        if chosen is None:
            return None
        ordered.append(pending.pop(chosen))
    return ordered


def _cut_crossing_word(contour, punctures, direction):
    """Derive the exact signed cut-crossing word of a based loop.

    This proof path intentionally accepts only contours for which SymPy can
    solve every cut-line intersection as a finite exact set and prove every
    retained crossing transverse.  Ambiguous/tangent cases are left unknown.
    """
    if not contour.closed:
        return None
    t = contour.parameter
    # Give real/imag expansion a parameter carrying an explicit real
    # assumption even when the caller's Symbol omitted one.
    tr = sp.Dummy("_t_real", real=True)
    gamma = contour.gamma.xreplace({t: tr})
    a = (
        contour.lower.xreplace({t: tr})
        if isinstance(contour.lower, sp.Basic)
        else contour.lower
    )
    b = (
        contour.upper.xreplace({t: tr})
        if isinstance(contour.upper, sp.Basic)
        else contour.upper
    )
    dc = sp.conjugate(sp.sympify(direction))
    interval = sp.Interval(sp.Min(a, b), sp.Max(a, b))
    events = []
    for j, p in enumerate(punctures, start=1):
        along, transverse = _real_imag((gamma - p) * dc)
        roots = sp.solveset(transverse, tr, domain=interval)
        if not isinstance(roots, sp.FiniteSet):
            return None
        dtrans = sp.diff(transverse, tr)
        for root in roots:
            # Endpoints are excluded because a valid cut system misses the
            # common basepoint.  If symbolic equality is undecidable, refuse.
            if (
                certify_identity_equal(root, a) is True
                or certify_identity_equal(root, b) is True
            ):
                continue
            if not _interior(root, a, b):
                return None
            radial = _strict_sign(along.subs(tr, root))
            if radial == 0:
                # The contour hits the puncture, so it is not a loop in the
                # punctured plane at all.
                return None
            if radial == -1:
                continue
            if radial is None:
                return None
            crossing = _strict_sign(dtrans.subs(tr, root))
            if crossing in (None, 0):
                return None
            # Parameter reversal reverses crossing orientation.
            if _strict_sign(sp.simplify(b - a)) == -1:
                crossing = -crossing
            events.append((root, j, crossing))
    ordered = _ordered_events(events, a, b)
    if ordered is None:
        return None
    return tuple(sign * generator for _, generator, sign in ordered)


def _free_reduce(word):
    stack = []
    for letter in word:
        letter = int(letter)
        if letter == 0:
            raise ValueError("free-group letters must be nonzero")
        if stack and stack[-1] == -letter:
            stack.pop()
        else:
            stack.append(letter)
    return tuple(stack)


_CUT_DIRECTIONS = (
    sp.S.One,
    sp.I,
    -sp.S.One,
    -sp.I,
    1 + sp.I,
    -1 + sp.I,
    -1 - sp.I,
    1 - sp.I,
    2 + sp.I,
    1 + 2 * sp.I,
)


def homotopic(contour1, contour2, *, punctures=(), return_certificate=False):
    """Certify based-loop homotopy in ``C`` minus finitely many punctures.

    The implementation uses a certified cut system and the standard free-group
    normal form of ``pi_1(C \\ P)``.  Thus generator order is retained and the
    result is genuinely nonabelian for two or more punctures.  ``None`` means
    that exact crossing geometry could not be certified for the supplied
    symbolic parameterizations; it is not a negative homotopy result.
    """
    if not (
        isinstance(contour1, ParametricContour)
        and isinstance(contour2, ParametricContour)
    ):
        return None
    if (
        certify_identity_equal(contour1.start, contour2.start) is not True
        or certify_identity_equal(contour1.end, contour2.end) is not True
    ):
        return False
    if not contour1.closed or not contour2.closed:
        return None
    ps = _distinct_punctures(punctures)
    # For zero or one puncture the free group is trivial/cyclic, so the exact
    # winding-number certificate is already complete and covers substantially
    # more symbolic contours than explicit cut-intersection solving.
    if len(ps) <= 1:
        if not ps:
            a = b = ()
        else:
            wa = winding_number(contour1, ps[0])
            wb = winding_number(contour2, ps[0])
            if wa is None or wb is None:
                return None
            if wa != wb:
                return False
            a = tuple(([1] * wa) if wa >= 0 else ([-1] * (-wa)))
            b = tuple(([1] * wb) if wb >= 0 else ([-1] * (-wb)))
        cert = HomotopyCertificate(contour1, contour2, ps, None, a, b, _free_reduce(a))
        return cert if return_certificate else True
    base = contour1.start
    for direction in _CUT_DIRECTIONS:
        if not _cut_system_valid(ps, base, direction):
            continue
        a = _cut_crossing_word(contour1, ps, direction)
        b = _cut_crossing_word(contour2, ps, direction)
        if a is None or b is None:
            continue
        ar = _free_reduce(a)
        br = _free_reduce(b)
        if ar != br:
            return False
        cert = HomotopyCertificate(contour1, contour2, ps, direction, a, b, ar)
        return cert if return_certificate else True
    return None


@dataclass(frozen=True)
class ContinuationCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    contour: ParametricContour
    continued: sp.Expr
    substitutions: tuple

    def verify(self):
        return (
            analytic_continue(
                self.expression, self.variable, self.contour, return_certificate=False
            )
            == self.continued
        )


def _divisor_winding(expr, variable, contour):
    """Return the (possibly rational) winding of an algebraic product/quotient.

    The result is the change of a continuous logarithm divided by ``2*pi*i``.
    Rational powers therefore retain their valuation
    as simple zeros.  ``None`` means that this exact structural proof does not
    cover the expression.
    """
    e = sp.sympify(expr)
    z = sp.sympify(variable)
    if not e.has(z):
        return sp.S.Zero if e.is_zero is False else None
    if e.is_Mul:
        parts = [_divisor_winding(arg, z, contour) for arg in e.args]
        if any(part is None for part in parts):
            return None
        return sp.simplify(sum(parts, sp.S.Zero))
    if e.is_Pow and e.exp.is_Rational is True:
        base_winding = _divisor_winding(e.base, z, contour)
        if base_winding is None:
            return None
        return sp.simplify(e.exp * base_winding)
    num, den = sp.fraction(sp.cancel(e))
    total = sp.S.Zero
    for poly, sign in ((num, 1), (den, -1)):
        try:
            P = sp.Poly(poly, z)
        except sp.PolynomialError:
            return None
        roots = sp.roots(P)
        if sum(roots.values()) != P.degree():
            return None
        for root, multiplicity in roots.items():
            winding = winding_number(contour, root)
            if winding is None:
                return None
            total += sign * multiplicity * winding
    return sp.simplify(total)


def analytic_continue(expr, var, contour, *, return_certificate=False):
    """Continue explicit logarithms and rational powers along a closed contour."""
    e = sp.sympify(expr)
    z = sp.sympify(var)
    if not contour.closed:
        return None
    repl = {}
    evidence = []
    atoms = sorted(e.atoms(sp.log), key=sp.default_sort_key)
    for a in atoms:
        arg = a.args[0]
        delta = _divisor_winding(arg, z, contour)
        if delta is None:
            return None
        repl[a] = a + 2 * sp.pi * sp.I * delta
        evidence.append((a, delta))
    # Principal rational powers acquire exp(2*pi*i*q*w) around zeros/poles of the base.
    for a in sorted(e.atoms(sp.Pow), key=sp.default_sort_key):
        if a.exp.is_Rational is not True or a.exp.is_integer is True:
            continue
        base = a.base
        num, den = sp.fraction(sp.cancel(base))
        delta = 0
        for poly, sgn in ((num, 1), (den, -1)):
            try:
                P = sp.Poly(poly, z)
            except sp.PolynomialError:
                return None
            roots = sp.roots(P)
            if sum(roots.values()) != P.degree():
                return None
            for root, mult in roots.items():
                w = winding_number(contour, root)
                if w is None:
                    return None
                delta += sgn * mult * w
        repl[a] = a * sp.exp(2 * sp.pi * sp.I * a.exp * delta)
        evidence.append((a, delta))
    out = sp.expand(e.xreplace(repl))
    c = ContinuationCertificate(e, z, contour, out, tuple(evidence))
    return c if return_certificate else out


def monodromy(expr, var, contour, *, return_certificate=False):
    return analytic_continue(expr, var, contour, return_certificate=return_certificate)


@dataclass(frozen=True)
class SingularityClassificationCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    point: sp.Expr
    kind: str
    order: sp.Expr | None

    def verify(self):
        fresh = classify_singularity(
            self.expression, self.variable, self.point, return_certificate=True
        )
        return (
            fresh is not None and fresh.kind == self.kind and fresh.order == self.order
        )


def classify_singularity(expr, var, point, *, return_certificate=False):
    """Classify an isolated point using exact local order and branch structure."""
    from .singularities.records import order_at

    e = sp.sympify(expr)
    z = sp.sympify(var)
    p = sp.sympify(point)
    # log(g(z)) has logarithmic monodromy when g has a nonzero local divisor
    # order, whether that order is a zero or a pole.
    for a in e.atoms(sp.log):
        ao = order_at(a.args[0], z, p)
        if ao is not None and ao != 0:
            c = SingularityClassificationCertificate(e, z, p, "branch_point", None)
            return c if return_certificate else c.kind
    # For g(z)^q, q rational nonintegral, branching is controlled by q times
    # the local divisor order.  sqrt(z**2), for example, has integral local
    # exponent and therefore no local monodromy at zero.
    for a in e.atoms(sp.Pow):
        if a.exp.is_Rational is True and a.exp.is_integer is not True:
            bo = order_at(a.base, z, p)
            if bo is not None and bo != 0:
                local_exp = sp.simplify(a.exp * bo)
                if local_exp.is_integer is not True:
                    c = SingularityClassificationCertificate(
                        e, z, p, "branch_point", None
                    )
                    return c if return_certificate else c.kind
    o = order_at(e, z, p)
    if o is not None and o.is_integer is True:
        kind = (
            "pole"
            if o.is_negative is True
            else ("zero" if o.is_positive is True else "regular")
        )
        c = SingularityClassificationCertificate(e, z, p, kind, o)
        return c if return_certificate else kind
    # exp/sin/cos/... of a pole has an essential singularity.
    for a in e.atoms(sp.exp, sp.sin, sp.cos, sp.sinh, sp.cosh):
        ao = order_at(a.args[0], z, p)
        if ao is not None and ao.is_negative is True:
            c = SingularityClassificationCertificate(e, z, p, "essential", None)
            return c if return_certificate else c.kind
    return None


@dataclass(frozen=True)
class ContourDeformationCertificate:
    expression: sp.Expr
    variable: sp.Symbol
    first: ParametricContour
    second: ParametricContour
    punctures: tuple

    def verify(self):
        return (
            deform_contour(
                self.expression,
                self.variable,
                self.first,
                self.second,
                punctures=self.punctures,
                return_certificate=False,
            )
            is True
        )


def deform_contour(
    expr, var, first, second, *, punctures=None, return_certificate=False
):
    """Certify contour deformation for rational meromorphic integrands."""
    e = sp.cancel(sp.sympify(expr))
    z = sp.sympify(var)
    _num, den = sp.fraction(e)
    try:
        pd = sp.Poly(den, z)
        roots = sp.roots(pd)
        if sum(roots.values()) != pd.degree():
            return None
        singular = tuple(roots)
    except sp.PolynomialError:
        return None
    # User-supplied punctures are additional forbidden points; they must never
    # replace poles forced by the integrand, or an omitted pole could make an
    # invalid deformation appear certified.
    extra = () if punctures is None else tuple(map(sp.sympify, punctures))
    ps = tuple(dict.fromkeys((*singular, *extra)))
    h = homotopic(first, second, punctures=ps)
    if h is not True:
        return h
    c = ContourDeformationCertificate(e, z, first, second, ps)
    return c if return_certificate else True
