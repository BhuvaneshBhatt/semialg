from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..engine import truth_bool
from ..exactness import certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS


@dataclass(frozen=True)
class AnalyticNeighborhoodCertificate:
    holomorphic: bool | None
    region: object
    method: str = "unknown"
    obstructions: tuple[object, ...] = ()
    evidence: tuple[object, ...] = ()

    @property
    def certified(self):
        return self.holomorphic is not None


_ENTIRE = (sp.exp, sp.sin, sp.cos, sp.sinh, sp.cosh, sp.erf, sp.erfc, sp.erfi)


def _all3(values):
    unknown = False
    for v in values:
        if v is False:
            return False
        if v is None:
            unknown = True
    return None if unknown else True


def _rect_contains(rect, z):
    a, b = rect
    z = sp.sympify(z)
    tests = (
        sp.Le(sp.re(a), sp.re(z)),
        sp.Le(sp.re(z), sp.re(b)),
        sp.Le(sp.im(a), sp.im(z)),
        sp.Le(sp.im(z), sp.im(b)),
    )
    vals = [truth_bool(t) for t in tests]
    if any(v is False for v in vals):
        return False
    if all(v is True for v in vals):
        return True
    return None


def _affine_data(g, var):
    try:
        p = sp.Poly(sp.expand(g), var)
        if p.degree() > 1:
            return None
        c = p.nth(1)
        d = p.nth(0)
        if c.has(var) or d.has(var):
            return None
        return c, d
    except SYMBOLIC_QUERY_ERRORS:
        return None


def _affine_image_bounds(g, var, rect):
    # exact rectangular enclosure for affine g(z)=c*z+d
    cd = _affine_data(g, var)
    if cd is None:
        return None
    c, d = cd
    try:
        a, b = rect
        corners = [
            sp.re(a) + sp.I * sp.im(a),
            sp.re(a) + sp.I * sp.im(b),
            sp.re(b) + sp.I * sp.im(a),
            sp.re(b) + sp.I * sp.im(b),
        ]
        vals = [sp.expand_complex(c * q + d) for q in corners]
        rs = [sp.re(v) for v in vals]
        ims = [sp.im(v) for v in vals]
        if not all(v.is_real is True for v in rs + ims):
            return None
        # default_sort_key is syntactic, not numerical ordering.  Min/Max keep
        # the bound exact and simplify when comparisons are decidable.
        return sp.Min(*rs), sp.Max(*rs), sp.Min(*ims), sp.Max(*ims)
    except SYMBOLIC_QUERY_ERRORS:
        return None


def _avoids_negative_real_axis(g, var, rect):
    b = _affine_image_bounds(g, var, rect)
    if b is None:
        return None
    r0, _r1, i0, i1 = b
    # definitely above/below real axis, or entirely right of 0
    if i0.is_positive is True or i1.is_negative is True or r0.is_positive is True:
        return True
    intersects_box = (
        r0.is_nonpositive is True
        and i0.is_nonpositive is True
        and i1.is_nonnegative is True
    )
    if not intersects_box:
        return None
    # Bounding-box intersection does not prove intersection of a rotated
    # parallelogram with the cut.  It *is* exact for real/pure-imaginary affine
    # maps, whose image of an axis-aligned rectangle is axis-aligned.
    cd = _affine_data(g, var)
    if cd is not None:
        c, _ = cd
        axis_aligned = c.is_real is True or certify_zero(sp.re(c)) is True
        if axis_aligned:
            return False
    return None


def _nonzero_on_rectangle(g, var, rect):
    if not g.has(var):
        z = sp.simplify(g)
        zero = certify_zero(z)
        if zero is True:
            return False
        if zero is False:
            return True
        return None
    try:
        sol = sp.solveset(sp.Eq(g, 0), var, domain=S.Complexes)
        if isinstance(sol, sp.FiniteSet):
            unknown = False
            for r in sol:
                q = _rect_contains(rect, r)
                if q is True:
                    return False
                if q is None:
                    unknown = True
            if not unknown:
                return True
    except SYMBOLIC_QUERY_ERRORS:
        pass
    b = _affine_image_bounds(g, var, rect)
    if b is not None:
        r0, r1, i0, i1 = b
        if (
            r0.is_positive is True
            or r1.is_negative is True
            or i0.is_positive is True
            or i1.is_negative is True
        ):
            return True
    return None


def holomorphic_on_rectangle(expr, var, rectangle):
    """Conservatively certify principal-branch holomorphy on a closed rectangle.

    Returns ``holomorphic=None`` when the current symbolic tests cannot certify
    either outcome.  No floating sampling is used for a positive certificate.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    rect = tuple(map(sp.sympify, rectangle))
    obs = []
    ev = []

    def rec(e):
        if not e.has(var) or e == var:
            return True
        if e.func in _ENTIRE:
            return _all3(rec(a) for a in e.args)
        if e.func in (sp.Add, sp.Mul):
            return _all3(rec(a) for a in e.args)
        if e.func is sp.Pow:
            b, p = e.args
            sub = _all3((rec(b), rec(p)))
            if sub is False:
                return False
            if p.is_Integer:
                if p.is_nonnegative:
                    return sub
                nz = _nonzero_on_rectangle(b, var, rect)
                ev.append(("nonzero", b, nz))
                if nz is False:
                    obs.append(sp.Eq(b, 0))
                    return False
                return _all3((sub, nz))
            av = _avoids_negative_real_axis(b, var, rect)
            ev.append(("principal_power_cut", b, av))
            if av is False:
                obs.append(("negative_real_axis", b))
                return False
            return _all3((sub, av))
        if e.func is sp.log:
            sub = rec(e.args[0])
            if sub is False:
                return False
            av = _avoids_negative_real_axis(e.args[0], var, rect)
            ev.append(("log_cut", e.args[0], av))
            if av is False:
                obs.append(("negative_real_axis", e.args[0]))
                return False
            return _all3((sub, av))
        if e.func in (sp.tan,):
            sub = rec(e.args[0])
            if sub is False:
                return False
            nz = _nonzero_on_rectangle(sp.cos(e.args[0]), var, rect)
            ev.append(("tan_poles", nz))
            if nz is False:
                obs.append(sp.Eq(sp.cos(e.args[0]), 0))
                return False
            return _all3((sub, nz))
        if e.func in (sp.gamma,):
            sub = rec(e.args[0])
            if sub is False:
                return False
            # finite exact preimages of nonpositive integers need a dedicated
            # pole-set proof; do not turn uncertainty into non-holomorphy.
            return None
        return None

    ans = rec(expr)
    return AnalyticNeighborhoodCertificate(
        ans, rect, "structural_rectangle", tuple(obs), tuple(ev)
    )


def analytic_neighborhood(expr, var, interval, radius=None):
    """Certify a complex rectangular neighborhood of a finite real interval."""
    interval = sp.sympify(interval)
    if (
        not isinstance(interval, sp.Interval)
        or interval.start in (-sp.oo, sp.oo)
        or interval.end in (-sp.oo, sp.oo)
    ):
        return AnalyticNeighborhoodCertificate(None, interval, "unsupported_interval")
    if radius is None:
        radius = sp.Rational(1, 100)
    r = sp.sympify(radius)
    if r.is_positive is not True:
        return AnalyticNeighborhoodCertificate(None, interval, "invalid_radius")
    rect = (interval.start - sp.I * r, interval.end + sp.I * r)
    return holomorphic_on_rectangle(expr, var, rect)
