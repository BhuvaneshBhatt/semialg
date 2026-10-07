from __future__ import annotations

import contextlib
import itertools
from dataclasses import dataclass

import mpmath as mp
import sympy as sp
from sympy import S

from ..engine import truth_bool
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from .analytic import AnalyticNeighborhoodCertificate, holomorphic_on_rectangle


@dataclass(frozen=True)
class ComplexRectangle:
    """Closed axis-aligned complex rectangle with exact real bounds."""

    re_min: sp.Expr
    re_max: sp.Expr
    im_min: sp.Expr
    im_max: sp.Expr

    @property
    def lower_left(self):
        return self.re_min + sp.I * self.im_min

    @property
    def upper_right(self):
        return self.re_max + sp.I * self.im_max

    @property
    def center(self):
        return (self.re_min + self.re_max) / 2 + sp.I * (self.im_min + self.im_max) / 2

    @property
    def width(self):
        return sp.simplify(self.re_max - self.re_min)

    @property
    def height(self):
        return sp.simplify(self.im_max - self.im_min)

    def as_pair(self):
        return self.lower_left, self.upper_right


@dataclass(frozen=True)
class ArgumentPrincipleCertificate:
    rectangle: ComplexRectangle
    winding_number: int | None
    zero_count: int | None
    pole_count: int | None = 0
    certified: bool = False
    method: str = "unknown"
    boundary_segments: int = 0
    precision: int | None = None
    analytic_certificate: AnalyticNeighborhoodCertificate | None = None
    evidence: tuple[object, ...] = ()


@dataclass(frozen=True)
class CertifiedComplexRootRectangle:
    rectangle: ComplexRectangle
    root_count: int
    certificate: ArgumentPrincipleCertificate


@dataclass(frozen=True)
class ComplexRootIsolationResult:
    roots: tuple[sp.Expr, ...] = ()
    rectangles: tuple[CertifiedComplexRootRectangle, ...] = ()
    complete: bool = False
    total_count: int | None = None
    method: str = "unknown"
    evidence: tuple[object, ...] = ()
    certificates: tuple[ArgumentPrincipleCertificate, ...] = ()

    @property
    def certified(self):
        return self.complete


def _normalize_rectangle(rectangle):
    if isinstance(rectangle, ComplexRectangle):
        rect = rectangle
    else:
        try:
            a, b = map(sp.sympify, rectangle)
        except (TypeError, ValueError):
            return None
        bounds = tuple(map(sp.simplify, (sp.re(a), sp.re(b), sp.im(a), sp.im(b))))
        if not all(q.is_real is True for q in bounds):
            return None
        rect = ComplexRectangle(
            sp.Min(bounds[0], bounds[1]),
            sp.Max(bounds[0], bounds[1]),
            sp.Min(bounds[2], bounds[3]),
            sp.Max(bounds[2], bounds[3]),
        )
    if truth_bool(rect.re_min < rect.re_max) is not True:
        return None
    if truth_bool(rect.im_min < rect.im_max) is not True:
        return None
    return rect


def _contains(rect: ComplexRectangle, root):
    root = sp.sympify(root)
    checks = (
        rect.re_min <= sp.re(root),
        sp.re(root) <= rect.re_max,
        rect.im_min <= sp.im(root),
        sp.im(root) <= rect.im_max,
    )
    vals = tuple(truth_bool(q) for q in checks)
    if all(v is True for v in vals):
        return True
    if any(v is False for v in vals):
        return False
    return None


def _exact_roots(expr, var, rect):
    with contextlib.suppress(*SYMBOLIC_QUERY_ERRORS):
        solved = sp.solveset(sp.Eq(expr, 0), var, domain=S.Complexes)
        if isinstance(solved, sp.FiniteSet):
            roots = []
            for root in solved:
                inside = _contains(rect, root)
                if inside is None:
                    return None
                if inside:
                    roots.append(root)
            return tuple(sorted(set(roots), key=sp.default_sort_key))
    return None


def _iv_modules():
    names = (
        "sin",
        "cos",
        "tan",
        "asin",
        "acos",
        "atan",
        "sinh",
        "cosh",
        "tanh",
        "asinh",
        "acosh",
        "atanh",
        "exp",
        "log",
        "sqrt",
        "erf",
        "erfc",
    )
    return {name: getattr(mp.iv, name) for name in names if hasattr(mp.iv, name)}


def _interval_function(expr, var):
    try:
        return sp.lambdify(var, expr, modules=[_iv_modules(), "mpmath"])
    except (TypeError, ValueError, NameError):
        return None


def _mp_number(value, dps):
    return mp.mpf(str(sp.N(value, dps + 10)))


def _segment_box(z0, z1, dps):
    r0, r1 = sorted((_mp_number(sp.re(z0), dps), _mp_number(sp.re(z1), dps)))
    i0, i1 = sorted((_mp_number(sp.im(z0), dps), _mp_number(sp.im(z1), dps)))
    return mp.iv.mpc([r0, r1], [i0, i1])


def _image_excludes_zero(value):
    try:
        return not (0 in value.real and 0 in value.imag)
    except (TypeError, ValueError, AttributeError):
        return False


def _point_value(function, z, dps):
    box = _segment_box(z, z, dps)
    value = function(box)
    if not hasattr(value, "real"):
        value = mp.iv.mpc(value)
    # The midpoint of an interval point evaluation is a valid representative of
    # its convex enclosure.  It is used only after the containing boundary image
    # has been certified zero-free.
    re = (mp.mpf(value.real.a) + mp.mpf(value.real.b)) / 2
    im = (mp.mpf(value.imag.a) + mp.mpf(value.imag.b)) / 2
    return mp.mpc(re, im)


def _certify_segment(function, z0, z1, *, dps, depth, max_depth, vertices):
    try:
        image = function(_segment_box(z0, z1, dps))
        if not hasattr(image, "real"):
            image = mp.iv.mpc(image)
    except (ArithmeticError, TypeError, ValueError, ZeroDivisionError):
        return False
    if _image_excludes_zero(image):
        vertices.append(z1)
        return True
    if depth >= max_depth:
        return False
    mid = sp.simplify((z0 + z1) / 2)
    return _certify_segment(
        function,
        z0,
        mid,
        dps=dps,
        depth=depth + 1,
        max_depth=max_depth,
        vertices=vertices,
    ) and _certify_segment(
        function,
        mid,
        z1,
        dps=dps,
        depth=depth + 1,
        max_depth=max_depth,
        vertices=vertices,
    )


def _boundary_vertices(rect, function, *, dps, max_depth):
    a = rect.re_min
    b = rect.re_max
    c = rect.im_min
    d = rect.im_max
    corners = (
        a + sp.I * c,
        b + sp.I * c,
        b + sp.I * d,
        a + sp.I * d,
        a + sp.I * c,
    )
    vertices = [corners[0]]
    for z0, z1 in itertools.pairwise(corners):
        if not _certify_segment(
            function,
            z0,
            z1,
            dps=dps,
            depth=0,
            max_depth=max_depth,
            vertices=vertices,
        ):
            return None
    return tuple(vertices)


def _polygon_winding(values):
    total = mp.mpf("0")
    for a, b in itertools.pairwise(values):
        if a == 0 or b == 0:
            return None
        cross = a.real * b.imag - a.imag * b.real
        dot = a.real * b.real + a.imag * b.imag
        total += mp.atan2(cross, dot)
    turns = total / (2 * mp.pi)
    winding = int(mp.nint(turns))
    if abs(turns - winding) > mp.mpf("1e-20"):
        return None
    return winding


def _numerical_winding(expr, var, rect, *, dps=60, samples=256):
    try:
        f = sp.lambdify(var, expr, "mpmath")
        a = complex(float(sp.N(rect.re_min)), float(sp.N(rect.im_min)))
        b = complex(float(sp.N(rect.re_max)), float(sp.N(rect.im_min)))
        c = complex(float(sp.N(rect.re_max)), float(sp.N(rect.im_max)))
        d = complex(float(sp.N(rect.re_min)), float(sp.N(rect.im_max)))
    except (TypeError, ValueError, NameError):
        return None
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        points = []
        per_edge = max(8, samples // 4)
        for start, stop in ((a, b), (b, c), (c, d), (d, a)):
            for k in range(per_edge):
                t = k / per_edge
                points.append(start + t * (stop - start))
        points.append(a)
        values = []
        for z in points:
            try:
                values.append(complex(f(z)))
            except (ArithmeticError, TypeError, ValueError, ZeroDivisionError):
                return None
        return _polygon_winding(values)
    finally:
        mp.mp.dps = old


def argument_principle_count(
    expr,
    var,
    rectangle,
    *,
    precision=60,
    max_boundary_depth=14,
):
    """Count zeros in a rectangle by a validated argument-principle contour.

    A positive certificate requires a symbolic proof of holomorphy on the closed
    rectangle and interval-arithmetic enclosures of every boundary segment that
    exclude zero.  If validation fails, a numerical winding estimate may be
    returned as evidence but ``certified`` remains false.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    rect = _normalize_rectangle(rectangle)
    if rect is None:
        return ArgumentPrincipleCertificate(
            ComplexRectangle(S.Zero, S.One, S.Zero, S.One),
            None,
            None,
            certified=False,
            method="invalid_rectangle",
        )
    analytic = holomorphic_on_rectangle(expr, var, rect.as_pair())
    function = _interval_function(expr, var)
    if analytic.holomorphic is True and function is not None:
        old_iv_dps = mp.iv.dps
        old_mp_dps = mp.mp.dps
        mp.iv.dps = int(precision)
        mp.mp.dps = int(precision)
        try:
            vertices = _boundary_vertices(
                rect,
                function,
                dps=int(precision),
                max_depth=max_boundary_depth,
            )
            if vertices is not None:
                try:
                    values = tuple(
                        _point_value(function, z, int(precision)) for z in vertices
                    )
                    winding = _polygon_winding(values)
                except (ArithmeticError, TypeError, ValueError, ZeroDivisionError):
                    winding = None
                if winding is not None and winding >= 0:
                    return ArgumentPrincipleCertificate(
                        rect,
                        winding,
                        winding,
                        0,
                        True,
                        "validated_interval_argument_principle",
                        len(vertices) - 1,
                        int(precision),
                        analytic,
                        (("boundary_zero_free", True),),
                    )
        finally:
            mp.iv.dps = old_iv_dps
            mp.mp.dps = old_mp_dps
    estimate = _numerical_winding(expr, var, rect, dps=int(precision))
    return ArgumentPrincipleCertificate(
        rect,
        estimate,
        None,
        0 if analytic.holomorphic is True else None,
        False,
        "numerical_argument_principle_estimate",
        0,
        int(precision),
        analytic,
        (("winding_estimate", estimate),),
    )


def _children(rect, rx, ry):
    x = sp.simplify(rect.re_min + rx * (rect.re_max - rect.re_min))
    y = sp.simplify(rect.im_min + ry * (rect.im_max - rect.im_min))
    return (
        ComplexRectangle(rect.re_min, x, rect.im_min, y),
        ComplexRectangle(x, rect.re_max, rect.im_min, y),
        ComplexRectangle(x, rect.re_max, y, rect.im_max),
        ComplexRectangle(rect.re_min, x, y, rect.im_max),
    )


def isolate_complex_roots(
    expr,
    var,
    rectangle,
    *,
    precision=60,
    max_boundary_depth=14,
    max_isolation_depth=10,
):
    """Return certified complex root rectangles inside ``rectangle``.

    Rectangles with count one isolate a unique simple zero.  A terminal rectangle
    with count greater than one is a certified cluster counted with multiplicity.
    Completeness means that the certified rectangle counts add to the certified
    argument-principle count of the original rectangle.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    rect = _normalize_rectangle(rectangle)
    if rect is None:
        return ComplexRootIsolationResult(method="invalid_rectangle")

    exact = _exact_roots(expr, var, rect)
    if exact is not None:
        return ComplexRootIsolationResult(
            exact,
            complete=True,
            total_count=len(exact),
            method="exact_complex_roots",
            evidence=(("rectangle", rect),),
        )

    cache = {}

    def count(r):
        key = (r.re_min, r.re_max, r.im_min, r.im_max)
        if key not in cache:
            cache[key] = argument_principle_count(
                expr,
                var,
                r,
                precision=precision,
                max_boundary_depth=max_boundary_depth,
            )
        return cache[key]

    parent = count(rect)
    if not parent.certified or parent.zero_count is None:
        return ComplexRootIsolationResult(
            complete=False,
            total_count=None,
            method="uncertified_complex_isolation",
            evidence=(parent,),
            certificates=(parent,),
        )
    if parent.zero_count == 0:
        return ComplexRootIsolationResult(
            complete=True,
            total_count=0,
            method="validated_argument_principle",
            evidence=(parent,),
            certificates=(parent,),
        )

    leaves = []
    certs = [parent]
    split_ratios = (
        (sp.Rational(1, 2), sp.Rational(1, 2)),
        (sp.Rational(2, 5), sp.Rational(2, 5)),
        (sp.Rational(3, 5), sp.Rational(3, 5)),
        (sp.Rational(2, 5), sp.Rational(3, 5)),
        (sp.Rational(3, 5), sp.Rational(2, 5)),
    )

    def rec(r, cert, depth):
        n = cert.zero_count
        if n is None or n == 0:
            return True
        if n == 1:
            leaves.append(CertifiedComplexRootRectangle(r, 1, cert))
            return True
        if depth >= max_isolation_depth:
            leaves.append(CertifiedComplexRootRectangle(r, n, cert))
            return True
        for rx, ry in split_ratios:
            child_rects = _children(r, rx, ry)
            child_certs = tuple(count(q) for q in child_rects)
            if not all(c.certified and c.zero_count is not None for c in child_certs):
                continue
            if sum(c.zero_count for c in child_certs) != n:
                continue
            certs.extend(child_certs)
            return all(
                rec(q, c, depth + 1)
                for q, c in zip(child_rects, child_certs, strict=True)
            )
        leaves.append(CertifiedComplexRootRectangle(r, n, cert))
        return True

    complete = rec(rect, parent, 0)
    return ComplexRootIsolationResult(
        rectangles=tuple(leaves),
        complete=complete,
        total_count=parent.zero_count,
        method="validated_complex_rectangle_isolation",
        evidence=(parent,),
        certificates=tuple(certs),
    )


def complex_zeros(
    expr,
    var,
    *,
    rectangle,
    return_result=False,
    precision=60,
    max_boundary_depth=14,
    max_isolation_depth=10,
):
    """Return zeros in a closed complex rectangle, with certification when needed."""
    result = isolate_complex_roots(
        expr,
        var,
        rectangle,
        precision=precision,
        max_boundary_depth=max_boundary_depth,
        max_isolation_depth=max_isolation_depth,
    )
    if return_result:
        return result
    if result.complete and result.roots:
        return sp.FiniteSet(*result.roots)
    if result.complete and result.total_count == 0:
        return S.EmptySet
    return result
