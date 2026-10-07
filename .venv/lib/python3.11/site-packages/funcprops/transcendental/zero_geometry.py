from __future__ import annotations

import itertools
from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..exactness import certify_equal, certify_identity_equal, certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS


class IndexedZero(sp.Expr):
    """Exact symbolic indexed zero of a named special-function family.

    The object is intentionally unevaluated: it carries exact identity/order
    information without pretending the zero has an elementary closed form.
    """

    is_real = True

    def __new__(cls, family, index, order=S.NaN, derivative_order=0, branch="real"):
        fam = sp.Symbol(str(family)) if not isinstance(family, sp.Symbol) else family
        idx = sp.sympify(index)
        ordv = sp.sympify(order)
        d = sp.Integer(derivative_order)
        br = sp.Symbol(str(branch))
        return sp.Expr.__new__(cls, fam, idx, ordv, d, br)

    @property
    def family(self):
        return str(self.args[0])

    @property
    def index(self):
        return self.args[1]

    @property
    def order(self):
        return None if self.args[2] is S.NaN else self.args[2]

    @property
    def derivative_order(self):
        return int(self.args[3])

    @property
    def branch(self):
        return str(self.args[4])

    def _sympystr(self, printer):
        fam = self.family
        k = printer.doprint(self.index)
        d = self.derivative_order
        if fam in ("BesselJ", "BesselY"):
            nu = printer.doprint(self.args[2])
            base = "j" if fam == "BesselJ" else "y"
            return f"{base}_{nu},{k}" + ("'" * d)
        if fam in ("AiryAi", "AiryBi"):
            base = "a" if fam == "AiryAi" else "b"
            return f"{base}_{k}" + ("'" * d)
        return f"Zero[{fam},{k},d={d}]"


@dataclass(frozen=True)
class ZeroOrdering:
    left: object
    right: object
    relation: str = "<"
    guard: sp.Expr = S.true
    source: str = ""


@dataclass(frozen=True)
class ZeroInterlacing:
    first_family: str
    second_family: str
    pattern: tuple[ZeroOrdering, ...]
    guard: sp.Expr = S.true
    source: str = ""
    complete: bool = True


@dataclass(frozen=True)
class ZeroGeometryCertificate:
    family: str
    support: sp.Set
    count: object
    guard: sp.Expr = S.true
    simple: bool | None = None
    indexed_family: tuple[IndexedZero, ...] = ()
    orderings: tuple[ZeroOrdering, ...] = ()
    interlacings: tuple[ZeroInterlacing, ...] = ()
    derivative_order: int = 0
    complete_on_support: bool = True
    source: str = ""
    notes: tuple[object, ...] = ()
    real_exhaustive: bool = False


@dataclass(frozen=True)
class IntervalPropertyCell:
    interval: object
    value: object
    guard: sp.Expr = S.true
    proof: tuple[object, ...] = ()


@dataclass(frozen=True)
class PropertyIntervalDecomposition:
    derivative_order: int
    domain: object
    cells: tuple[IntervalPropertyCell, ...] = ()
    boundary_geometry: tuple[ZeroGeometryCertificate, ...] = ()
    alternating: bool | None = None
    first_value: object = None
    complete: bool = False
    source: str = ""


def _is_named(expr, name):
    return getattr(getattr(expr, "func", None), "__name__", "") == name


def _bessel_geometry(expr, var, derivative_order=0):
    f = expr.func
    if f not in (sp.besselj, sp.bessely) or len(expr.args) != 2 or expr.args[1] != var:
        return None
    nu = expr.args[0]
    fam = "BesselJ" if f is sp.besselj else "BesselY"
    support = sp.Interval.open(0, sp.oo)
    # DLMF 10.21(i): infinitely many positive real zeros for real order;
    # derivative simplicity has the sharper order guards stated there.
    if derivative_order == 0:
        guard = sp.Contains(nu, S.Reals)
        simple = True
    elif fam == "BesselJ":
        guard = sp.And(sp.Contains(nu, S.Reals), nu >= -1)
        simple = True
    else:
        guard = sp.And(sp.Contains(nu, S.Reals), nu >= -sp.Rational(1, 2))
        simple = True
    k = sp.Symbol("k", integer=True, positive=True)
    z = IndexedZero(fam, k, nu, derivative_order, "positive")
    orders = []
    inter = []
    if derivative_order == 0:
        # Adjacent orders interlace for nu >= 0, DLMF 10.21.E2.
        a = IndexedZero(fam, k, nu, 0, "positive")
        b = IndexedZero(fam, k, nu + 1, 0, "positive")
        c = IndexedZero(fam, k + 1, nu, 0, "positive")
        rel = (
            ZeroOrdering(a, b, "<", nu >= 0, "DLMF 10.21.E2"),
            ZeroOrdering(b, c, "<", nu >= 0, "DLMF 10.21.E2"),
        )
        inter.append(
            ZeroInterlacing(
                fam + f"({nu})", fam + f"({nu + 1})", rel, nu >= 0, "DLMF 10.21.E2"
            )
        )
        # Function/derivative zeros alternate; E3 gives the combined chain for nu>=0.
        dp = IndexedZero(fam, k, nu, 1, "positive")
        dn = IndexedZero(fam, k + 1, nu, 1, "positive")
        rel2 = (
            ZeroOrdering(dp, a, "<", nu >= 0, "DLMF 10.21.E3"),
            ZeroOrdering(a, dn, "<", nu >= 0, "DLMF 10.21.E3"),
        )
        inter.append(ZeroInterlacing(fam + "'", fam, rel2, nu >= 0, "DLMF 10.21.E3"))
    src = "DLMF 10.21(i), 10.21.E2-E3"
    return ZeroGeometryCertificate(
        fam,
        support,
        sp.oo,
        guard,
        simple,
        (z,),
        tuple(orders),
        tuple(inter),
        derivative_order,
        True,
        src,
    )


def _airy_geometry(expr, var, derivative_order=0):
    mapping = {
        sp.airyai: ("AiryAi", 0),
        sp.airyaiprime: ("AiryAi", 1),
        sp.airybi: ("AiryBi", 0),
        sp.airybiprime: ("AiryBi", 1),
    }
    if expr.func not in mapping or len(expr.args) != 1 or expr.args[0] != var:
        return None
    fam, builtin_d = mapping[expr.func]
    d = builtin_d + derivative_order
    if d not in (0, 1):
        return None
    support = sp.Interval.open(-sp.oo, 0)
    k = sp.Symbol("k", integer=True, positive=True)
    z = IndexedZero(fam, k, S.NaN, d, "negative")
    if fam == "AiryAi":
        # a'_{k+1} < a_k < a'_k
        rel = (
            ZeroOrdering(
                IndexedZero(fam, k + 1, S.NaN, 1, "negative"),
                IndexedZero(fam, k, S.NaN, 0, "negative"),
                "<",
                S.true,
                "DLMF 9.9 + Sturm separation",
            ),
            ZeroOrdering(
                IndexedZero(fam, k, S.NaN, 0, "negative"),
                IndexedZero(fam, k, S.NaN, 1, "negative"),
                "<",
                S.true,
                "DLMF 9.9 + Sturm separation",
            ),
        )
    else:
        # b_{k+1} < b'_k < b_k
        rel = (
            ZeroOrdering(
                IndexedZero(fam, k + 1, S.NaN, 0, "negative"),
                IndexedZero(fam, k, S.NaN, 1, "negative"),
                "<",
                S.true,
                "DLMF 9.9 + Sturm separation",
            ),
            ZeroOrdering(
                IndexedZero(fam, k, S.NaN, 1, "negative"),
                IndexedZero(fam, k, S.NaN, 0, "negative"),
                "<",
                S.true,
                "DLMF 9.9 + Sturm separation",
            ),
        )
    inter = ZeroInterlacing(
        fam, fam + "'", rel, S.true, "DLMF 9.9(i), E1-E4 + Sturm separation"
    )
    return ZeroGeometryCertificate(
        fam,
        support,
        sp.oo,
        S.true,
        True,
        (z,),
        (),
        (inter,),
        d,
        True,
        "DLMF 9.9(i), E1-E4",
        real_exhaustive=True,
    )


def _orthogonal_geometry(expr, var):
    if not expr.args or expr.args[-1] != var:
        return None
    f = expr.func
    n = expr.args[0]
    n_guard = sp.And(sp.Contains(n, S.Integers), n >= 1)
    if f is sp.jacobi and len(expr.args) == 4:
        a, b = expr.args[1:3]
        support = sp.Interval.open(-1, 1)
        guard = sp.And(n_guard, a > -1, b > -1)
        fam = "Jacobi"
        src = "DLMF 18.16.E1 + 18.2(vi)"
    elif f is sp.gegenbauer and len(expr.args) == 3:
        lam = expr.args[1]
        support = sp.Interval.open(-1, 1)
        guard = sp.And(n_guard, lam > -sp.Rational(1, 2), sp.Ne(lam, 0))
        fam = "Gegenbauer"
        src = "DLMF 18.16(iii) + 18.2(vi)"
    elif f is sp.assoc_laguerre and len(expr.args) == 3:
        a = expr.args[1]
        support = sp.Interval.open(0, sp.oo)
        guard = sp.And(n_guard, a > -1)
        fam = "Laguerre"
        src = "DLMF 18.16.E9 + 18.2(vi)"
    elif f is sp.laguerre and len(expr.args) == 2:
        support = sp.Interval.open(0, sp.oo)
        guard = n_guard
        fam = "Laguerre"
        src = "DLMF 18.16.E9 + 18.2(vi)"
    elif f in tuple(
        q
        for q in (getattr(sp, "hermite", None), getattr(sp, "hermite_prob", None))
        if q is not None
    ):
        support = S.Reals
        guard = n_guard
        fam = "Hermite"
        src = "DLMF 18.16(v) + 18.2(vi)"
    elif f in (sp.chebyshevt, sp.chebyshevu) and len(expr.args) == 2:
        support = sp.Interval.open(-1, 1)
        guard = n_guard
        fam = "ChebyshevT" if f is sp.chebyshevt else "ChebyshevU"
        src = "DLMF 18.16(iii) + 18.2(vi)"
    else:
        return None
    k = sp.Symbol("k", integer=True, positive=True)
    z = IndexedZero(fam, k, n, 0, "support")
    return ZeroGeometryCertificate(
        fam,
        support,
        n,
        guard,
        True,
        (z,),
        (),
        (),
        0,
        True,
        src,
        ("consecutive degrees interlace",),
        True,
    )


def _legendre_geometry(expr, var):
    # SymPy's assoc_legendre is Ferrers P on (-1,1) in the integer-degree use case.
    if expr.func is sp.assoc_legendre and len(expr.args) == 3 and expr.args[2] == var:
        nu, mu = expr.args[:2]
        support = sp.Interval.open(-1, 1)
        guard = sp.And(
            sp.Contains(nu, S.Naturals0), sp.Contains(mu, S.Integers), sp.Abs(mu) <= nu
        )
        count = sp.Max(nu - sp.Abs(mu), 0)
        return ZeroGeometryCertificate(
            "FerrersP",
            support,
            count,
            guard,
            True,
            (),
            (),
            (),
            0,
            True,
            "DLMF 14.16(ii)",
        )
    # DLMF-native Q may appear through catalogue/user-defined SymPy Function.
    if _is_named(expr, "LegendreQ") or _is_named(expr, "FerrersQ"):
        args = expr.args
        if len(args) >= 3 and args[-1] == var:
            nu, mu = args[0], args[1]
        elif len(args) >= 2 and args[-1] == var:
            nu, mu = args[0], S.Zero
        else:
            return None
        support = sp.Interval.open(-1, 1)
        base = sp.Max(sp.ceiling(nu - sp.Abs(mu)), 0)
        # DLMF 14.16(ii): Q zeros interlace P; count is base+k, k in {-1,0,1,2}
        # subject to a parity condition. Preserve this as an exact finite choice set.
        count = sp.FiniteSet(base - 1, base, base + 1, base + 2)
        guard = sp.And(
            sp.Contains(nu, S.Reals), sp.Contains(mu, S.Reals), nu >= -sp.Rational(1, 2)
        )
        IndexedZero(
            "FerrersP", sp.Symbol("k", integer=True, positive=True), nu, 0, "(-1,1)"
        )
        q = IndexedZero(
            "FerrersQ", sp.Symbol("k", integer=True, positive=True), nu, 0, "(-1,1)"
        )
        inter = ZeroInterlacing(
            "FerrersP", "FerrersQ", (), guard, "DLMF 14.16(ii): Q zeros interlace P"
        )
        return ZeroGeometryCertificate(
            "FerrersQ",
            support,
            count,
            guard,
            True,
            (q,),
            (),
            (inter,),
            0,
            True,
            "DLMF 14.16(ii)",
            ("count parity depends on cos(pi*nu), cos(pi*mu)",),
        )
    return None


def _zero_geometry_impl(expr, var, *, derivative_order=0, domain=S.Reals):
    """Return exact zero-geometry certificates known for ``expr`` or its derivative.

    This layer does not numerically locate zeros.  It records indexed families,
    counts, simplicity, support, and interlacing/order relations.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    # Strip nonzero variable-independent scalar factors.
    if expr.func is sp.Mul:
        c, dep = expr.as_independent(var, as_Add=False)
        if c != 0 and dep != expr:
            return _zero_geometry_impl(
                dep, var, derivative_order=derivative_order, domain=domain
            )
    # Elementary complete indexed families share the same ZeroGeometry contract.
    if derivative_order == 0:
        try:
            from ..infinite_root_families import structured_zero_geometry

            sg = structured_zero_geometry(expr, var)
        except (TypeError, ValueError, NotImplementedError):
            sg = None
        if sg is not None:
            return (sg,)
    # For requested derivative geometry, prefer the base family's theorem.
    for fn in (_bessel_geometry, _airy_geometry):
        g = fn(expr, var, derivative_order)
        if g is not None:
            return (g,)
    if derivative_order == 0:
        for fn in (_orthogonal_geometry, _legendre_geometry):
            g = fn(expr, var)
            if g is not None:
                return (g,)
    # Generic derivative closure: many orthogonal families differentiate into
    # another member of the same hierarchy.  Recurse only on the explicit
    # derivative expression, never numerically.
    if derivative_order:
        try:
            from ..proof_cache import cached_derivative

            d = cached_derivative(expr, var, order=derivative_order, factored=True)
        except SYMBOLIC_QUERY_ERRORS:
            d = None
        if d is not None and d != expr:
            return _zero_geometry_impl(d, var, derivative_order=0, domain=domain)
    return ()


def zero_geometry(expr, var, *, derivative_order=0, domain=S.Reals):
    """Cached exact zero-geometry proof primitive."""
    from ..proof_cache import property_proof_context

    return property_proof_context(expr, var, domain=domain).zero_geometry(
        derivative_order=derivative_order
    )


def _interval_subset_exact(a, b):
    """Prove ``a`` is a subset of interval ``b`` from exact endpoint facts."""
    if not isinstance(a, sp.Interval) or not isinstance(b, sp.Interval):
        return None
    # left endpoint
    if b.start is -sp.oo:
        left = True
    elif a.start is -sp.oo:
        left = False
    else:
        if (
            b.start == 0
            and isinstance(a.start, IndexedZero)
            and a.start.branch == "positive"
        ):
            left = True
        elif (
            b.start == 0
            and isinstance(a.start, IndexedZero)
            and a.start.branch == "negative"
        ):
            left = False
        else:
            d = sp.simplify(a.start - b.start)
            if d.is_positive is True:
                left = True
            elif certify_equal(a.start, b.start) is True:
                left = not (b.left_open and not a.left_open)
            elif d.is_negative is True:
                left = False
            else:
                left = None
    # right endpoint
    if b.end is sp.oo:
        right = True
    elif a.end is sp.oo:
        right = False
    else:
        if b.end == 0 and isinstance(a.end, IndexedZero) and a.end.branch == "negative":
            right = True
        elif (
            b.end == 0 and isinstance(a.end, IndexedZero) and a.end.branch == "positive"
        ):
            right = False
        else:
            d = sp.simplify(b.end - a.end)
            if d.is_positive is True:
                right = True
            elif certify_equal(b.end, a.end) is True:
                right = not (b.right_open and not a.right_open)
            elif d.is_negative is True:
                right = False
            else:
                right = None
    if left is False or right is False:
        return False
    if left is True and right is True:
        return True
    return None


def geometry_complete_on_domain(g, domain):
    """Whether a zero-geometry certificate exhausts zeros on ``domain``."""
    if not g.complete_on_support:
        return False
    try:
        if domain == g.support or domain.is_subset(g.support) is True:
            return True
    except SYMBOLIC_QUERY_ERRORS:
        pass
    if _interval_subset_exact(domain, g.support) is True:
        return True
    if g.real_exhaustive:
        try:
            if domain is S.Reals or domain.is_subset(S.Reals) is True:
                return True
        except SYMBOLIC_QUERY_ERRORS:
            pass
    return False


def geometry_support_inside_domain(g, domain):
    """Whether the certified indexed support is contained in ``domain``."""
    try:
        if (
            domain is S.Reals
            or g.support == domain
            or g.support.is_subset(domain) is True
        ):
            return True
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return _interval_subset_exact(g.support, domain) is True


def property_interval_decomposition(
    expr,
    var,
    *,
    derivative_order=1,
    domain=S.Reals,
    positive_value=None,
    negative_value=None,
    derivative=None,
):
    """Build an exact interval-pattern certificate from derivative zero geometry.

    Infinite indexed boundaries are represented by ``boundary_geometry`` and an
    alternating flag.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    geoms = zero_geometry(expr, var, derivative_order=derivative_order, domain=domain)
    if geoms:
        simple = all(g.simple is True for g in geoms)

        def guard_true(g):
            if g.guard is S.true:
                return True
            try:
                return sp.simplify(g.guard) is S.true
            except SYMBOLIC_QUERY_ERRORS:
                return False

        complete = all(
            geometry_complete_on_domain(g, domain) and guard_true(g) for g in geoms
        )
        # A simple zero of the derivative changes its sign.  Thus the property
        # (monotonicity for order 1; curvature for order 2) alternates by cell.
        return PropertyIntervalDecomposition(
            derivative_order,
            domain,
            (),
            geoms,
            True if simple else None,
            None,
            complete,
            "certified derivative-zero geometry",
        )
    # Generic finite exact derivative roots: materialize concrete cells.
    if isinstance(domain, sp.Interval):
        try:
            from ..proof_cache import cached_derivative

            d = (
                sp.factor(derivative)
                if derivative is not None
                else cached_derivative(expr, var, order=derivative_order, factored=True)
            )
            sol = sp.solveset(sp.Eq(d, 0), var, domain=domain)
        except SYMBOLIC_QUERY_ERRORS:
            sol = None
        if sol is S.EmptySet or isinstance(sol, sp.FiniteSet):
            pts = []
            for q in () if sol is S.EmptySet else sol:
                placed = False
                for i, r in enumerate(pts):
                    delta = sp.simplify(q - r)
                    if delta.is_negative is True:
                        pts.insert(i, q)
                        placed = True
                        break
                    if delta.is_positive is True:
                        continue
                    if certify_equal(q, r) is True:
                        placed = True
                        break
                    pts = None
                    break
                if pts is None:
                    break
                if not placed:
                    pts.append(q)
            if pts is not None:
                # Critical points at interval endpoints do not split the open
                # interior into additional monotonic cells.  Keeping them here
                # creates empty cells such as Interval.open(a,a) and can turn a
                # complete monotonicity proof into an artificial failure.
                pts = [q for q in pts if q != domain.start and q != domain.end]
                simple = True
                try:
                    from ..proof_cache import cached_derivative

                    dn = cached_derivative(d, var, order=1)
                    for q in pts:
                        if certify_zero(dn.subs(var, q)) is not False:
                            simple = False
                            break
                except SYMBOLIC_QUERY_ERRORS:
                    simple = False
                bounds = [domain.start, *pts, domain.end]
                cells = []
                complete = True
                for a, b in itertools.pairwise(bounds):
                    if a is -sp.oo and b is sp.oo:
                        sample = S.Zero
                    elif a is -sp.oo:
                        sample = b - 1
                    elif b is sp.oo:
                        sample = a + 1
                    else:
                        sample = sp.simplify((a + b) / 2)
                    sg = sp.simplify(d.subs(var, sample))
                    if sg.is_positive is True:
                        val = positive_value
                    elif sg.is_negative is True:
                        val = negative_value
                    else:
                        complete = False
                        val = None
                    cells.append(
                        IntervalPropertyCell(
                            sp.Interval.open(a, b),
                            val,
                            S.true,
                            (("derivative_sample", sample, sg),),
                        )
                    )
                fg = ZeroGeometryCertificate(
                    "finite-derivative-zeros",
                    sp.FiniteSet(*pts),
                    len(pts),
                    S.true,
                    simple,
                    (),
                    (),
                    (),
                    derivative_order,
                    True,
                    "exact solveset + derivative multiplicity test",
                )
                return PropertyIntervalDecomposition(
                    derivative_order,
                    domain,
                    tuple(cells),
                    (fg,),
                    simple if pts else False,
                    None,
                    complete,
                    "finite exact derivative-root decomposition",
                )

    # Airy second derivative: y'' = x*y.  This supplies exact curvature
    # boundaries from the function zeros plus x=0.
    if (
        derivative_order == 2
        and expr.func in (sp.airyai, sp.airybi)
        and expr.args == (var,)
    ):
        base = zero_geometry(expr, var, domain=domain)
        if base:
            z0 = ZeroGeometryCertificate(
                "linear-factor-x",
                sp.FiniteSet(0),
                1,
                S.true,
                True,
                (IndexedZero("x", 1, S.NaN, 0, "zero"),),
                (),
                (),
                0,
                True,
                "Airy ODE y''=x y",
            )
            return PropertyIntervalDecomposition(
                2,
                domain,
                (),
                (*base, z0),
                True,
                positive_value,
                True,
                "Airy ODE + certified Airy zeros",
            )
    return PropertyIntervalDecomposition(
        derivative_order,
        domain,
        (),
        (),
        None,
        None,
        False,
        "insufficient zero geometry",
    )


def zero_ordering_bounds(zero):
    """Return exact symbolic neighbor inequalities known for an ``IndexedZero``."""
    if not isinstance(zero, IndexedZero):
        return ()
    k = zero.index
    d = zero.derivative_order
    fam = zero.family
    out = []
    if fam in ("BesselJ", "BesselY") and zero.order is not None:
        nu = zero.order
        if d == 0:
            # DLMF 10.21.E2 (adjacent order) and E3 (function/derivative chain).
            out.append(
                ZeroOrdering(
                    zero,
                    IndexedZero(fam, k, nu + 1, 0, "positive"),
                    "<",
                    nu >= 0,
                    "DLMF 10.21.E2",
                )
            )
            out.append(
                ZeroOrdering(
                    IndexedZero(fam, k, nu, 1, "positive"),
                    zero,
                    "<",
                    nu >= 0,
                    "DLMF 10.21.E3",
                )
            )
            out.append(
                ZeroOrdering(
                    zero,
                    IndexedZero(fam, k + 1, nu, 1, "positive"),
                    "<",
                    nu >= 0,
                    "DLMF 10.21.E3",
                )
            )
            out.append(
                ZeroOrdering(
                    IndexedZero(fam, k - 1, nu + 1, 0, "positive"),
                    zero,
                    "<",
                    sp.And(nu >= 0, k >= 2),
                    "DLMF 10.21.E2",
                )
            )
        elif d == 1:
            out.append(
                ZeroOrdering(
                    zero,
                    IndexedZero(fam, k, nu, 0, "positive"),
                    "<",
                    nu >= 0,
                    "DLMF 10.21.E3",
                )
            )
    elif fam == "AiryAi":
        if d == 0:
            out.extend(
                (
                    ZeroOrdering(
                        IndexedZero(fam, k + 1, S.NaN, 1, "negative"),
                        zero,
                        "<",
                        S.true,
                        "DLMF 9.9 + Sturm separation",
                    ),
                    ZeroOrdering(
                        zero,
                        IndexedZero(fam, k, S.NaN, 1, "negative"),
                        "<",
                        S.true,
                        "DLMF 9.9 + Sturm separation",
                    ),
                )
            )
        elif d == 1:
            out.extend(
                (
                    ZeroOrdering(
                        IndexedZero(fam, k, S.NaN, 0, "negative"),
                        zero,
                        "<",
                        S.true,
                        "DLMF 9.9 + Sturm separation",
                    ),
                    ZeroOrdering(
                        zero,
                        IndexedZero(fam, k - 1, S.NaN, 0, "negative"),
                        "<",
                        k >= 2,
                        "DLMF 9.9 + Sturm separation",
                    ),
                )
            )
    elif fam == "AiryBi":
        if d == 1:
            out.extend(
                (
                    ZeroOrdering(
                        IndexedZero(fam, k + 1, S.NaN, 0, "negative"),
                        zero,
                        "<",
                        S.true,
                        "DLMF 9.9 + Sturm separation",
                    ),
                    ZeroOrdering(
                        zero,
                        IndexedZero(fam, k, S.NaN, 0, "negative"),
                        "<",
                        S.true,
                        "DLMF 9.9 + Sturm separation",
                    ),
                )
            )
    return tuple(out)


@dataclass(frozen=True)
class CertifiedIntegerBracket:
    """Rigorous integer enclosure for an exact indexed-tail selector."""

    lower: object
    upper: object
    guard: sp.Expr = S.true
    source: str = ""
    proof: tuple[object, ...] = ()

    @property
    def exact(self):
        return certify_identity_equal(self.lower, self.upper) is True

    def as_interval(self):
        return sp.Interval(sp.sympify(self.lower), sp.sympify(self.upper))


def _besselj_zero_linear_bounds(order, derivative_order=0):
    """Return rigorous affine-in-k bounds for positive J zeros.

    For nu>=0, DLMF 10.21(iv) gives monotonicity in order and 10.21.E2
    gives adjacent-order interlacing.  Anchoring at J_{-1/2}(x), whose kth
    positive zero is (k-1/2)pi, yields

        (k-1/2)pi <= j_{nu,k} < (k+N-1/2)pi,
        N = ceiling(nu+1/2).

    DLMF 10.21.E3 then gives, for k>=2,
        j_{nu,k-1} < j'_{nu,k} < j_{nu,k}.
    """
    nu = sp.sympify(order)
    if nu.is_nonnegative is not True:
        return None
    k = sp.Symbol("_k", integer=True, positive=True)
    # Upper anchor: choose the first half-integer order 1/2+M >= nu.
    # Repeated DLMF 10.21.E2 interlacing gives
    # j_{1/2+M,k} < j_{1/2,k+M}=(k+M)pi.
    M = sp.Max(0, sp.ceiling(nu - sp.Rational(1, 2)))
    uj = (k + M) * sp.pi
    # Lower anchor is always J_-1/2; when nu>=1/2 we can strengthen it
    # to J_1/2, whose kth zero is k*pi.
    half_or_more = (nu - sp.Rational(1, 2)).is_nonnegative is True
    base_lj = k * sp.pi if half_or_more else (k - sp.Rational(1, 2)) * sp.pi
    if derivative_order == 0:
        return (
            k,
            base_lj,
            uj,
            sp.Ge(nu, 0),
            (
                "DLMF 10.21(iv) order monotonicity",
                "DLMF 10.21.E2 interlacing",
                "J_±1/2 elementary zeros",
            ),
        )
    if derivative_order == 1:
        prev_lj = (k - 1) * sp.pi if half_or_more else (k - sp.Rational(3, 2)) * sp.pi
        lj = sp.Piecewise((S.Zero, sp.Eq(k, 1)), (prev_lj, True))
        return (
            k,
            lj,
            uj,
            sp.Ge(nu, 0),
            (
                "DLMF 10.21.E3 derivative/function interlacing",
                "DLMF 10.21(iv) order monotonicity",
                "J_±1/2 elementary zeros",
            ),
        )
    return None


def indexed_tail_integer_bracket(
    zero_template, cutoff, *, direction="right", inclusive=True
):
    """Certify an integer bracket for the first included indexed zero.

    Currently formula-level rigorous for positive Bessel-J function and first
    derivative zero families with nonnegative real order.  Unsupported families
    return ``None``.
    """
    if not isinstance(zero_template, IndexedZero):
        raise TypeError("zero_template must be an IndexedZero")
    A = sp.sympify(cutoff)
    if direction != "right" or zero_template.branch != "positive":
        return None
    if zero_template.family != "BesselJ":
        return None
    data = _besselj_zero_linear_bounds(
        zero_template.args[2], zero_template.derivative_order
    )
    if data is None:
        return None
    _k, _L, _U, guard, proof = data
    nu = zero_template.args[2]
    M = sp.Max(0, sp.ceiling(nu - sp.Rational(1, 2)))
    # U(k)=(k+M)pi. If U(k)<A (or <=A for a strict tail), k is
    # certainly excluded. L(k)>=A (or >A) certifies inclusion.
    qlo = A / sp.pi - M
    lower = sp.Max(1, sp.ceiling(qlo)) if inclusive else sp.Max(1, sp.floor(qlo) + 1)
    half_or_more = (nu - sp.Rational(1, 2)).is_nonnegative is True
    if zero_template.derivative_order == 0:
        shift = S.Zero if half_or_more else sp.Rational(1, 2)
        qhi = A / sp.pi + shift
        upper = sp.Max(1, sp.ceiling(qhi) if inclusive else sp.floor(qhi) + 1)
    else:
        shift = S.One if half_or_more else sp.Rational(3, 2)
        qhi = A / sp.pi + shift
        upper = sp.Max(2, sp.ceiling(qhi) if inclusive else sp.floor(qhi) + 1)
    return CertifiedIntegerBracket(
        sp.simplify(lower),
        sp.simplify(upper),
        guard,
        "rigorous Bessel-J zero bounds derived from DLMF 10.21",
        proof,
    )


class IndexedTailIndex(sp.Expr):
    """Exact inverse index for a monotone indexed zero family.

    ``direction='right'`` means the least positive integer ``k`` whose zero is
    at or to the right of ``cutoff``.  ``direction='left'`` is the analogous
    first zero encountered moving left.  The object is exact and intentionally
    unevaluated when no closed integer formula is known.
    """

    is_integer = True
    is_positive = True

    def __new__(
        cls,
        family,
        cutoff,
        order=S.NaN,
        derivative_order=0,
        direction="right",
        inclusive=True,
        branch="real",
    ):
        fam = sp.Symbol(str(family))
        cutoff = sp.sympify(cutoff)
        order = sp.sympify(order)
        d = sp.Integer(derivative_order)
        direction = sp.Symbol(str(direction))
        inc = S.true if inclusive else S.false
        branch = sp.Symbol(str(branch))
        return sp.Expr.__new__(cls, fam, cutoff, order, d, direction, inc, branch)

    @property
    def family(self):
        return str(self.args[0])

    @property
    def cutoff(self):
        return self.args[1]

    @property
    def order(self):
        return None if self.args[2] is S.NaN else self.args[2]

    @property
    def derivative_order(self):
        return int(self.args[3])

    @property
    def direction(self):
        return str(self.args[4])

    @property
    def inclusive(self):
        return self.args[5] is S.true

    @property
    def branch(self):
        return str(self.args[6])

    def indexed_zero(self):
        return IndexedZero(
            self.family, self, self.args[2], self.derivative_order, self.branch
        )

    def _sympystr(self, printer):
        op = (
            ">="
            if self.direction == "right" and self.inclusive
            else ">"
            if self.direction == "right"
            else "<="
            if self.inclusive
            else "<"
        )
        return f"first_k[{self.family}: zero_k {op} {printer.doprint(self.cutoff)}]"


@dataclass(frozen=True)
class IndexedTailInversionCertificate:
    family: str
    cutoff: object
    direction: str
    inclusive: bool
    first_index: object
    zero_template: IndexedZero
    guard: sp.Expr = S.true
    complete: bool = True
    source: str = ""
    proof: tuple[object, ...] = ()
    integer_bracket: object = None

    @property
    def first_zero(self):
        return IndexedZero(
            self.family,
            self.first_index,
            self.zero_template.args[2],
            self.zero_template.derivative_order,
            self.zero_template.branch,
        )


def _same_zero_sequence(a, b):
    return (
        isinstance(a, IndexedZero)
        and isinstance(b, IndexedZero)
        and a.family == b.family
        and a.args[2] == b.args[2]
        and a.derivative_order == b.derivative_order
        and a.branch == b.branch
    )


def indexed_tail_inversion(zero_template, cutoff, *, direction="right", inclusive=True):
    """Invert an indexed monotone zero family exactly at a cutoff.

    The routine first applies exact same-family/interlacing identities.  For a
    general cutoff it returns an ``IndexedTailIndex`` whose semantics are the
    certified least included index; no floating root finder is used.
    """
    if not isinstance(zero_template, IndexedZero):
        raise TypeError("zero_template must be an IndexedZero")
    cutoff = sp.sympify(cutoff)
    fam = zero_template.family
    d = zero_template.derivative_order
    ordv = zero_template.args[2]
    branch = zero_template.branch
    if direction not in ("right", "left"):
        raise ValueError("direction must be 'right' or 'left'")

    # Positive Bessel zero sequences increase with k; negative Airy sequences
    # decrease with k (indices run away from the origin).
    increasing = branch == "positive"
    decreasing = branch == "negative"

    # Exact inversion at a member of the same indexed family.
    if isinstance(cutoff, IndexedZero) and _same_zero_sequence(zero_template, cutoff):
        m = cutoff.index
        if (increasing and direction == "right") or (
            decreasing and direction == "left"
        ):
            idx = m if inclusive else m + 1
        elif (increasing and direction == "left") or (
            decreasing and direction == "right"
        ):
            idx = m if inclusive else sp.Max(1, m - 1)
        else:
            idx = IndexedTailIndex(fam, cutoff, ordv, d, direction, inclusive, branch)
        return IndexedTailInversionCertificate(
            fam,
            cutoff,
            direction,
            inclusive,
            sp.simplify(idx),
            zero_template,
            S.true,
            True,
            "same indexed zero family",
            ("strict monotonicity of indexed zeros in k",),
        )

    # Bessel J/Y function-vs-derivative interlacing, DLMF 10.21.E3:
    # z'_{nu,m} < z_{nu,m} < z'_{nu,m+1} for the supported nonnegative orders.
    if (
        isinstance(cutoff, IndexedZero)
        and fam in ("BesselJ", "BesselY")
        and (
            cutoff.family == fam
            and cutoff.args[2] == ordv
            and d == 1
            and cutoff.derivative_order == 0
        )
    ):
        m = cutoff.index
        idx = m + 1 if direction == "right" else m
        guard = sp.Ge(ordv, 0) if ordv is not S.NaN else S.true
        return IndexedTailInversionCertificate(
            fam,
            cutoff,
            direction,
            inclusive,
            sp.simplify(idx),
            zero_template,
            guard,
            True,
            "DLMF 10.21.E3",
            ("derivative/function zeros interlace",),
        )

    # Airy function/derivative interlacing.  For Ai:
    # a'_{m+1} < a_m < a'_m.  Bi has b_{m+1} < b'_m < b_m.
    if (
        isinstance(cutoff, IndexedZero)
        and fam in ("AiryAi", "AiryBi")
        and cutoff.family == fam
    ):
        m = cutoff.index
        if d == 1 and cutoff.derivative_order == 0:
            if fam == "AiryAi":
                idx = m if direction == "right" else m + 1
            else:  # b_{m+1} < b'_m < b_m
                idx = m if direction == "left" else sp.Max(1, m - 1)
            return IndexedTailInversionCertificate(
                fam,
                cutoff,
                direction,
                inclusive,
                sp.simplify(idx),
                zero_template,
                S.true,
                True,
                "DLMF 9.9 + Sturm separation",
                ("Airy function/derivative zeros interlace",),
            )

    idx = IndexedTailIndex(fam, cutoff, ordv, d, direction, inclusive, branch)
    # The selector itself is an exact inverse of the certified monotone indexed
    # family. Rigorous zero bounds may further collapse it to a concrete integer.
    src = "exact symbolic inversion of certified indexed-zero ordering"
    bracket = indexed_tail_integer_bracket(
        zero_template, cutoff, direction=direction, inclusive=inclusive
    )
    if bracket is not None and bracket.exact:
        idx = sp.sympify(bracket.lower)
        src += " + certified integer zero-bound inversion"
    return IndexedTailInversionCertificate(
        fam,
        cutoff,
        direction,
        inclusive,
        idx,
        zero_template,
        S.true,
        True,
        src,
        (
            "no numerical root approximation",
            "least included positive integer index",
        ),
        bracket,
    )
