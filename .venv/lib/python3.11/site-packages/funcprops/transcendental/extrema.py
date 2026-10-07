from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..exactness import certify_equal, certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from .zero_geometry import (
    IndexedZero,
    ZeroGeometryCertificate,
    geometry_complete_on_domain,
    indexed_tail_inversion,
    zero_geometry,
)


class IndexedExtremumValue(sp.Expr):
    """Exact value of a function at an indexed critical point.

    This is an unevaluated proof object, not a numerical approximation.  The
    underlying exact SymPy value is available through ``.value``.
    """

    is_real = True

    def __new__(cls, family, index, point, value, kind="extremum"):
        fam = sp.Symbol(str(family))
        idx = sp.sympify(index)
        point = sp.sympify(point)
        value = sp.sympify(value)
        kind = sp.Symbol(str(kind))
        return sp.Expr.__new__(cls, fam, idx, point, value, kind)

    @property
    def family(self):
        return str(self.args[0])

    @property
    def index(self):
        return self.args[1]

    @property
    def point(self):
        return self.args[2]

    @property
    def value(self):
        return self.args[3]

    @property
    def kind(self):
        return str(self.args[4])

    def _sympystr(self, printer):
        return f"Extremum[{self.family},{printer.doprint(self.index)}]={printer.doprint(self.value)}"


class ExtremumEnvelope(sp.Expr):
    """Exact infimum/supremum over an indexed extremum family.

    ``direction`` is ``'inf'`` or ``'sup'``.  The object intentionally stays
    symbolic when no theorem proves which indexed member realizes the bound.
    """

    is_real = True

    def __new__(cls, direction, template, index, index_guard=S.true):
        d = sp.Symbol(str(direction))
        template = sp.sympify(template)
        index = sp.sympify(index)
        guard = sp.sympify(index_guard)
        return sp.Expr.__new__(cls, d, template, index, guard)

    @property
    def direction(self):
        return str(self.args[0])

    @property
    def template(self):
        return self.args[1]

    @property
    def index(self):
        return self.args[2]

    @property
    def index_guard(self):
        return self.args[3]

    def _sympystr(self, printer):
        op = "inf" if self.direction == "inf" else "sup"
        return f"{op}_{{{printer.doprint(self.index)}: {printer.doprint(self.index_guard)}}}({printer.doprint(self.template)})"


@dataclass(frozen=True)
class ExtremumFamilyCertificate:
    family: str
    domain: object
    critical_geometry: tuple[ZeroGeometryCertificate, ...]
    extrema: tuple[IndexedExtremumValue, ...] = ()
    alternating: bool | None = None
    guard: sp.Expr = S.true
    complete: bool = False
    source: str = ""
    notes: tuple[object, ...] = ()


@dataclass(frozen=True)
class ExtremumOrderingTheorem:
    """Exact theorem ordering successive indexed extrema.

    ``magnitude_direction='decreasing'`` means ``|E_{k+1}| < |E_k|``
    on ``index_guard``.  ``global_max_index``/``global_min_index`` identify
    a member when the alternating sign pattern plus magnitude ordering proves
    first-member dominance.  The proof field records the exact Lyapunov/Sonin
    quantity used by funcprops; it is not numerical evidence.
    """

    family: str
    magnitude_direction: str
    guard: sp.Expr = S.true
    index_guard: sp.Expr = S.true
    global_max_index: object = None
    global_min_index: object = None
    proof: object = None
    source: str = ""
    notes: tuple[object, ...] = ()


@dataclass(frozen=True)
class ExtremumRange:
    """Certified symbolic interval image from endpoint and critical-point data."""

    lower: object
    upper: object
    lower_attained: bool | None
    upper_attained: bool | None
    domain: object
    complete: bool
    candidates: tuple[object, ...] = ()
    evidence: tuple[object, ...] = ()
    source: str = ""

    def as_condition(self, y):
        y = sp.sympify(y)
        if not self.complete:
            return None
        lo = (
            y >= self.lower
            if self.lower_attained is True
            else y > self.lower
            if self.lower_attained is False
            else None
        )
        hi = (
            y <= self.upper
            if self.upper_attained is True
            else y < self.upper
            if self.upper_attained is False
            else None
        )
        if lo is None or hi is None:
            return None
        return sp.And(lo, hi)


def _domain_interval(domain, var):
    if isinstance(domain, sp.Interval):
        return domain
    try:
        from ..real_analysis.common import condition_to_set

        return condition_to_set(domain, var)
    except SYMBOLIC_QUERY_ERRORS:
        return None


def _endpoint_candidate(expr, var, endpoint, side, attained):
    from ..proof_cache import cached_endpoint_limit

    val = cached_endpoint_limit(expr, var, endpoint, direction=side, attained=attained)
    if val is None or getattr(val, "has", lambda *_: False)(var):
        return None
    return (sp.simplify(val), bool(attained))


def extremum_geometry(expr, var, *, domain=S.Reals, partition=None):
    """Return exact indexed critical-point/extremum families for ``expr``."""
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if partition is None:
        try:
            from ..real_analysis.monotonic_certificates import monotonic_partition

            partition = monotonic_partition(
                expr, var, domain=domain, return_result=True
            )
        except SYMBOLIC_QUERY_ERRORS:
            partition = None
    geoms = ()
    if partition is not None and getattr(partition, "indexed", None) is not None:
        geoms = tuple(partition.indexed.boundary_geometry)
    if not geoms:
        geoms = zero_geometry(expr, var, derivative_order=1, domain=domain)
    if not geoms:
        return ()
    out = []
    for g in geoms:
        extrema = []
        for z in g.indexed_family:
            val = sp.simplify(expr.xreplace({var: z}))
            kind = "extremum"
            # DLMF 9.9.E4/E5 fixes the Airy-Ai extrema sign parity.
            if g.family == "AiryAi" and z.derivative_order == 1:
                kind = "max_if_odd_min_if_even"
            elif g.family == "AiryBi" and z.derivative_order == 1:
                kind = "min_if_odd_max_if_even"
            extrema.append(IndexedExtremumValue(g.family, z.index, z, val, kind))
        out.append(
            ExtremumFamilyCertificate(
                g.family,
                domain,
                (g,),
                tuple(extrema),
                alternating=True if g.simple is True else None,
                guard=g.guard,
                complete=geometry_complete_on_domain(g, domain),
                source=g.source,
                notes=("critical points are derivative zeros",),
            )
        )
    return tuple(out)


def _contains_support(domain, support):
    """Three-valued exact support containment; only True is accepted."""
    try:
        r = domain.is_superset(support)
        if r is True or r is S.true:
            return True
        if r is False or r is S.false:
            return False
    except SYMBOLIC_QUERY_ERRORS:
        pass
    try:
        d = sp.simplify(support - domain)
        if d is S.EmptySet:
            return True
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return None


def extremum_ordering_theorems(expr, var):
    """Return exact ordering theorems for indexed extremum families.

    The Airy proof is derived from ``E=y'^2-x*y^2``: ``E'=-y^2`` and
    ``(E/(-x))'=y'^2/x^2>0`` for ``x<0``.  At a critical point this ratio is
    ``y^2``, hence amplitudes decrease as the Airy zero index increases away
    from zero.

    For Bessel J, with ``E=x^2*y'^2+(x^2-nu^2)*y^2``, Bessel's equation gives
    ``E'=2*x*y^2`` and
    ``(E/(x^2-nu^2))'=-2*x^3*y'^2/(x^2-nu^2)^2 < 0`` for ``x>nu``.  At
    critical points the ratio is ``y^2``.  DLMF 10.21.E3 places the positive
    J_nu derivative zeros at or beyond nu (strictly beyond for the nonzero
    critical points), yielding decreasing successive extremum magnitudes.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if (
        expr.func in (sp.airyai, sp.airybi)
        and len(expr.args) == 1
        and expr.args[0] == var
    ):
        fam = "AiryAi" if expr.func is sp.airyai else "AiryBi"
        # DLMF 9.9.E4 fixes sign parity at derivative zeros.
        if fam == "AiryAi":
            gmax, gmin = sp.Integer(1), sp.Integer(2)
        else:
            gmax, gmin = sp.Integer(2), sp.Integer(1)
        return (
            ExtremumOrderingTheorem(
                fam,
                "decreasing",
                S.true,
                sp.Symbol("k", integer=True, positive=True) >= 1,
                gmax,
                gmin,
                "Airy energy: E=y'^2-x*y^2; (E/(-x))'=y'^2/x^2>0 on x<0; at y'=0, E/(-x)=y^2",
                "derived from Airy ODE + DLMF 9.9.E4; DLMF 9.8/9.9",
                ("index k increases away from zero, so |extremum_{k+1}|<|extremum_k|",),
            ),
        )
    if (
        expr.func in (sp.besselj, sp.bessely)
        and len(expr.args) == 2
        and expr.args[1] == var
    ):
        nu = expr.args[0]
        k = sp.Symbol("k", integer=True, positive=True)
        fam = "BesselJ" if expr.func is sp.besselj else "BesselY"
        # For nu>0: first critical point is a positive maximum, second a
        # negative minimum.  nu=0 uses the DLMF convention j'_{0,1}=0, so the
        # first *interior* critical point is index 2 and is a minimum.
        return (
            ExtremumOrderingTheorem(
                fam,
                "decreasing",
                sp.And(sp.Contains(nu, S.Reals), nu >= 0),
                k >= 1,
                sp.Piecewise((3, sp.Eq(nu, 0)), (1, True)),
                sp.Piecewise((2, sp.Eq(nu, 0)), (2, True)),
                "Bessel energy: E=x^2*y'^2+(x^2-nu^2)*y^2; (E/(x^2-nu^2))'=-2*x^3*y'^2/(x^2-nu^2)^2<0 for x>nu",
                "derived from Bessel ODE + DLMF 10.21.E1,E3; J bound additionally uses DLMF 10.14.E1",
                (
                    "for nu=0, j'_{0,1}=0 is an endpoint critical point, not a positive interior extremum",
                ),
            ),
        )
    return ()


def _indexed_value(expr, var, family, index, order=S.NaN, derivative_order=1):
    z = IndexedZero(
        family,
        sp.sympify(index),
        order,
        derivative_order,
        "positive" if family.startswith("Bessel") else "negative",
    )
    return sp.simplify(expr.xreplace({var: z})), z


def _first_odd_index(k):
    k = sp.sympify(k)
    return sp.simplify(k + sp.Mod(k + 1, 2))


def _first_even_index(k):
    k = sp.sympify(k)
    return sp.simplify(k + sp.Mod(k, 2))


def _tail_inversion_for_extrema(expr, var, domain):
    """Exact first critical-index certificate for an oscillatory semi-infinite tail."""
    if not isinstance(domain, sp.Interval):
        return None
    certs = extremum_geometry(expr, var, domain=domain)
    if not certs or not certs[0].extrema:
        return None
    z = certs[0].extrema[0].point
    if not isinstance(z, IndexedZero):
        return None
    if z.branch == "positive" and domain.end is sp.oo and domain.start is not -sp.oo:
        return indexed_tail_inversion(
            z, domain.start, direction="right", inclusive=not domain.left_open
        )
    if z.branch == "negative" and domain.start is -sp.oo and domain.end is not sp.oo:
        return indexed_tail_inversion(
            z, domain.end, direction="left", inclusive=not domain.right_open
        )
    return None


def _simplify_indexed_tail(expr, var, domain, certs, end_candidates):
    """Collapse an indexed envelope on a truncated oscillatory tail.

    The first included critical index is exact even when represented by an
    ``IndexedTailIndex``.  Alternating sign plus decreasing extremum magnitude
    then reduces the infinite family to the first included odd/even extrema.
    """
    inv = _tail_inversion_for_extrema(expr, var, domain)
    if inv is None:
        return None
    evidence = list(certs) + list(extremum_ordering_theorems(expr, var)) + [inv]
    k0 = sp.sympify(inv.first_index)

    if expr.func in (sp.besselj, sp.bessely):
        nu = expr.args[0]
        if nu.is_nonnegative is not True:
            return None
        # For J_0, DLMF counts j'_{0,1}=0.  Any strictly positive tail starts
        # with an interior critical index >=2.
        if (
            expr.func is sp.besselj
            and certify_zero(nu) is True
            and domain.start.is_positive is True
        ):
            k0 = sp.Max(2, k0)
        kmax = _first_odd_index(k0)
        kmin = _first_even_index(k0)
        fam = "BesselJ" if expr.func is sp.besselj else "BesselY"
        vmax, zmax = _indexed_value(expr, var, fam, kmax, nu)
        vmin, zmin = _indexed_value(expr, var, fam, kmin, nu)
        candidates = [
            ("critical", zmin, vmin, True),
            ("critical", zmax, vmax, True),
            *list(end_candidates),
        ]
        finite_ends = [
            (v, att) for side, pt, v, att in end_candidates if pt not in (-sp.oo, sp.oo)
        ]
        lo = sp.Min(vmin, *(v for v, _ in finite_ends)) if finite_ends else vmin
        hi = sp.Max(vmax, *(v for v, _ in finite_ends)) if finite_ends else vmax
        # A cutoff can lie just beyond a larger excluded extremum, so f(A) may
        # dominate the first same-sign included extremum.  On a closed tail it
        # is attained; on an open tail exact attainment depends on that symbolic
        # comparison and remains unknown.
        if domain.left_open and finite_ends:
            la = ua = None
        else:
            la = ua = True
        return ExtremumRange(
            lo,
            hi,
            la,
            ua,
            domain,
            True,
            tuple(candidates),
            tuple(evidence),
            "indexed-tail inversion + Bessel extremum amplitude ordering",
        )

    if expr.func in (sp.airyai, sp.airybi):
        # Oscillatory Airy tails run to -infinity; indices increase leftward.
        if domain.start is not -sp.oo:
            return None
        fam = "AiryAi" if expr.func is sp.airyai else "AiryBi"
        if fam == "AiryAi":
            kmax = _first_odd_index(k0)
            kmin = _first_even_index(k0)
        else:
            kmin = _first_odd_index(k0)
            kmax = _first_even_index(k0)
        vmax, zmax = _indexed_value(expr, var, fam, kmax)
        vmin, zmin = _indexed_value(expr, var, fam, kmin)
        # The -infinity Airy limit is 0 and is strictly between the first
        # included positive/negative extrema, so it cannot determine a bound.
        vals = [(v, att) for _, pt, v, att in end_candidates if pt is not -sp.oo]
        lo = sp.Min(vmin, *(v for v, _ in vals)) if vals else vmin
        hi = sp.Max(vmax, *(v for v, _ in vals)) if vals else vmax
        # If the finite endpoint belongs to the interval, every candidate in
        # the Min/Max is attained.  For an open cutoff, whether its limiting
        # value dominates is intentionally left unresolved.
        if domain.right_open:
            la = ua = None
        else:
            la = ua = True
        candidates = [
            ("critical", zmin, vmin, True),
            ("critical", zmax, vmax, True),
            *list(end_candidates),
        ]
        return ExtremumRange(
            lo,
            hi,
            la,
            ua,
            domain,
            True,
            tuple(candidates),
            tuple(evidence),
            "indexed-tail inversion + Airy extremum amplitude ordering",
        )
    return None


def _simplify_indexed_extrema(expr, var, domain, certs, end_candidates):
    """Apply exact first-extremum dominance theorems to an indexed envelope."""
    thms = extremum_ordering_theorems(expr, var)
    if not thms:
        return None
    vals = [(v, att, side) for side, _, v, att in end_candidates]
    t = thms[0]
    evidence = [*list(certs), t]

    # Airy: theorem covers every negative critical point.  It is safe to choose
    # the first positive and first negative extrema only when the requested
    # domain contains the entire negative critical support.
    if expr.func is sp.airyai:
        support = sp.Interval.open(-sp.oo, 0)
        if _contains_support(domain, support) is not True:
            return None
        fam = "AiryAi"
        vmax, zmax = _indexed_value(expr, var, fam, t.global_max_index)
        vmin, zmin = _indexed_value(expr, var, fam, t.global_min_index)
        lo = vmin
        hi = vmax
        la = True
        ua = True
        candidates = [
            ("critical", zmin, vmin, True),
            ("critical", zmax, vmax, True),
            *list(end_candidates),
        ]
        # Endpoint/limit candidates can dominate a signed family member.
        for v, att, _side in vals:
            d = sp.simplify(v - hi)
            if d.is_positive is True or certify_zero(d) is True:
                hi = v
                ua = att
            d = sp.simplify(v - lo)
            if d.is_negative is True or certify_zero(d) is True:
                lo = v
                la = att
        # Known Airy right-half-line behavior makes these comparisons exact
        # even when SymPy cannot order a symbolic indexed value against Ai(0).
        if expr.func is sp.airyai:
            # Ai is positive for x>=0 and decreases to 0; first negative-side
            # maximum/minimum dominate all x>=0.
            if domain.end is sp.oo or (
                domain.end.is_nonnegative is True if domain.end.is_finite else False
            ):
                pass
        elif expr.func is sp.airybi:
            # On x>=0 Bi is positive and increasing.  If +infinity is included,
            # the upper range is unbounded; on an endpoint 0 the endpoint Bi(0)
            # exceeds all later positive extrema on the negative side only if
            # its comparison can be proved below.
            if domain.end is sp.oo:
                hi = sp.oo
                ua = False
        return ExtremumRange(
            lo,
            hi,
            la,
            ua,
            domain,
            True,
            tuple(candidates),
            tuple(evidence),
            "indexed extrema simplified by exact Airy amplitude ordering",
        )

    # Bessel J: currently simplify the canonical positive semiaxis, where all
    # positive derivative zeros are present.  More restricted tails keep the
    # generic envelope because selecting the first included index is a separate
    # integer-inversion problem.
    if expr.func is sp.besselj:
        support = sp.Interval.open(0, sp.oo)
        if _contains_support(domain, support) is not True:
            return None
        nu = expr.args[0]
        if certify_zero(nu) is True:
            vmin, zmin = _indexed_value(expr, var, "BesselJ", 2, nu)
            lo = vmin
            la = True
            # J_0(0)=1 and DLMF 10.14.E1 gives |J_0(x)|<=1.  On (0,oo) the
            # upper endpoint is a nonattained limit; on [0,oo) it is attained.
            hi = S.One
            ua = bool(domain.contains(0) is S.true or domain.contains(0) is True)
            cand = [("critical", zmin, vmin, True), *list(end_candidates)]
            return ExtremumRange(
                lo,
                hi,
                la,
                ua,
                domain,
                True,
                tuple(cand),
                tuple(evidence),
                "Bessel J0 first-interior-extremum theorem",
            )
        if nu.is_positive is True:
            vmax, zmax = _indexed_value(expr, var, "BesselJ", 1, nu)
            vmin, zmin = _indexed_value(expr, var, "BesselJ", 2, nu)
            return ExtremumRange(
                vmin,
                vmax,
                True,
                True,
                domain,
                True,
                (
                    ("critical", zmin, vmin, True),
                    ("critical", zmax, vmax, True),
                    *tuple(end_candidates),
                ),
                tuple(evidence),
                "Bessel J_nu first-extrema dominance for nu>0",
            )
    if expr.func is sp.bessely:
        support = sp.Interval.open(0, sp.oo)
        if _contains_support(domain, support) is not True:
            return None
        nu = expr.args[0]
        if nu.is_nonnegative is True:
            vmax, zmax = _indexed_value(expr, var, "BesselY", 1, nu)
            # Y_nu -> -infinity at 0+ for nu>=0 (Y_0 logarithmically).
            return ExtremumRange(
                -sp.oo,
                vmax,
                False,
                True,
                domain,
                True,
                (("critical", zmax, vmax, True), *tuple(end_candidates)),
                tuple(evidence),
                "Bessel Y_nu first-maximum dominance for nu>=0",
            )
    # Bi on the whole real line: all real critical amplitudes decrease away
    # from zero; b'_1 is the first (negative) minimum, while Bi has no
    # nonnegative zeros and grows to +infinity.
    if expr.func is sp.airybi and domain.start is -sp.oo and domain.end is sp.oo:
        fam = "AiryBi"
        vmin, zmin = _indexed_value(expr, var, fam, 1)
        return ExtremumRange(
            vmin,
            sp.oo,
            True,
            False,
            domain,
            True,
            (("critical", zmin, vmin, True), *tuple(end_candidates)),
            tuple(evidence),
            "Airy Bi global minimum from first-extremum dominance",
        )
    return None


def simplify_extremum_range(expr, var, result):
    """Simplify an ``ExtremumRange`` using exact indexed-extremum theorems."""
    if not isinstance(result, ExtremumRange) or not result.complete:
        return result
    if not isinstance(result.domain, sp.Interval):
        return result
    certs = tuple(
        e for e in result.evidence if isinstance(e, ExtremumFamilyCertificate)
    )
    # Existing result evidence generally stores certificates directly, while
    # Accept the compact tuple form returned by extremum_geometry.
    if not certs:
        certs = extremum_geometry(expr, var, domain=result.domain)
    ends = tuple(
        c
        for c in result.candidates
        if isinstance(c, tuple) and c and c[0] in ("left", "right")
    )
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    # Preserve the strongest closed-form canonical-semiaxis theorems first.
    # Tail inversion is for genuinely truncated supports.
    simp = _simplify_indexed_extrema(expr, var, result.domain, certs, ends)
    if simp is not None:
        return simp
    tail = _simplify_indexed_tail(expr, var, result.domain, certs, ends)
    return tail if tail is not None else result


def _indexed_envelopes(domain, certs):
    """Build exact inf/sup envelopes over complete indexed critical families."""
    if not certs or not all(c.complete for c in certs):
        return None
    templates = []
    for c in certs:
        if not c.extrema:
            return None
        e = c.extrema[0]
        k = e.index
        # Filter indexed critical points to the requested interval exactly.
        point = e.point
        guard = sp.And(k >= 1, c.guard)
        if isinstance(domain, sp.Interval):
            if domain.start is not -sp.oo:
                guard = sp.And(
                    guard,
                    point > domain.start if domain.left_open else point >= domain.start,
                )
            if domain.end is not sp.oo:
                guard = sp.And(
                    guard,
                    point < domain.end if domain.right_open else point <= domain.end,
                )
        templates.append((e.value, k, guard, e))
    # Multiple derivative families are rare; preserve them as Min/Max of exact envelopes.
    infs = [ExtremumEnvelope("inf", v, k, g) for v, k, g, _ in templates]
    sups = [ExtremumEnvelope("sup", v, k, g) for v, k, g, _ in templates]
    return (
        sp.Min(*infs) if len(infs) > 1 else infs[0],
        sp.Max(*sups) if len(sups) > 1 else sups[0],
        tuple(t[3] for t in templates),
    )


def symbolic_extremum_range(expr, var, domain=S.Reals, *, partition=None):
    """Compute an exact range envelope from critical-point geometry.

    Closed bounded intervals with finitely many exact critical points reduce to
    ordinary exact ``Min``/``Max`` endpoints.  Infinite oscillatory critical
    families remain represented by ``ExtremumEnvelope`` objects.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    ds = _domain_interval(domain, var)
    if not isinstance(ds, sp.Interval):
        return None
    if partition is None:
        try:
            from ..real_analysis.monotonic_certificates import monotonic_partition

            partition = monotonic_partition(expr, var, domain=ds, return_result=True)
        except SYMBOLIC_QUERY_ERRORS:
            partition = None

    # Critical points are owned by the shared monotonic partition.  Sign-cell
    # completeness is not required for extrema: the separate flag certifies
    # that the derivative-root set itself is exhaustive.
    pts = None
    if (
        partition is not None
        and getattr(partition, "indexed", None) is None
        and getattr(partition, "critical_points_complete", False)
    ):
        pts = tuple(getattr(partition, "critical_points", ()))
    if (
        pts is not None
        and ds.start not in (-sp.oo, sp.oo)
        and ds.end not in (-sp.oo, sp.oo)
    ):
        candidates = []
        values = []
        lc = _endpoint_candidate(expr, var, ds.start, "+", not ds.left_open)
        rc = _endpoint_candidate(expr, var, ds.end, "-", not ds.right_open)
        if lc is None or rc is None:
            return None
        candidates.extend((("left", ds.start, *lc), ("right", ds.end, *rc)))
        for q in pts:
            try:
                v = sp.simplify(expr.subs(var, q))
            except SYMBOLIC_QUERY_ERRORS:
                return None
            values.append((v, True))
            candidates.append(("critical", q, v, True))
        values = [lc, rc, *values]
        lo = sp.Min(*(v for v, _ in values))
        hi = sp.Max(*(v for v, _ in values))

        def attainment(bound, want_min):
            # First try to prove that an attained candidate realizes the bound.
            for v, att in values:
                if not att:
                    continue
                if certify_equal(bound, v) is True:
                    return True
            # Prove non-attainment only when the bound is an open-endpoint limit
            # and every attained candidate lies strictly inside it.
            open_hits = []
            for v, att in values:
                if att:
                    continue
                if certify_equal(bound, v) is True:
                    open_hits.append(v)
            if open_hits:
                ok = True
                for v, att in values:
                    if not att:
                        continue
                    d = sp.simplify(v - bound)
                    if want_min and d.is_positive is not True:
                        ok = False
                        break
                    if not want_min and d.is_negative is not True:
                        ok = False
                        break
                if ok:
                    return False
            return None

        la = attainment(lo, True)
        ua = attainment(hi, False)
        return ExtremumRange(
            lo,
            hi,
            la,
            ua,
            ds,
            True,
            tuple(candidates),
            (("critical_points", pts),),
            "finite exact critical-point range",
        )

    certs = extremum_geometry(expr, var, domain=ds, partition=partition)
    env = _indexed_envelopes(ds, certs)
    if env is None:
        return None
    crit_lo, crit_hi, templates = env
    end_candidates = []
    left = (
        _endpoint_candidate(expr, var, ds.start, "+", not ds.left_open)
        if ds.start is not -sp.oo
        else _endpoint_candidate(expr, var, -sp.oo, "+", False)
    )
    right = (
        _endpoint_candidate(expr, var, ds.end, "-", not ds.right_open)
        if ds.end is not sp.oo
        else _endpoint_candidate(expr, var, sp.oo, "-", False)
    )
    if left is not None:
        end_candidates.append(("left", ds.start, *left))
    if right is not None:
        end_candidates.append(("right", ds.end, *right))
    vals = [v for _, _, v, _ in end_candidates]
    lo = sp.Min(crit_lo, *vals) if vals else crit_lo
    hi = sp.Max(crit_hi, *vals) if vals else crit_hi
    # Infinite envelopes may be attained by a member or only approached; unless
    # a stronger theorem resolves this, preserve the endpoint-attainment state.
    raw = ExtremumRange(
        lo,
        hi,
        None,
        None,
        ds,
        True,
        tuple(end_candidates) + templates,
        tuple(certs),
        "indexed derivative-zero extremum envelope",
    )
    return simplify_extremum_range(expr, var, raw)
