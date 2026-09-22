"""Exact winding-fill construction from arbitrary planar polygonal paths."""

from __future__ import annotations

from collections import defaultdict
from functools import cmp_to_key

import sympy as sp

from ._zero_testing import certified_equal
from .boundary_topology import PolygonalComponent, PolygonalSet
from .exact_arithmetic import compare_exact_reals
from .standard_regions import Polygon


def _eq(a, b):
    return certified_equal(sp.expand(a - b), 0) is True


def _sgn(a):
    return compare_exact_reals(sp.expand(a), sp.Integer(0))


def _cross(a, b):
    return sp.expand(a[0] * b[1] - a[1] * b[0])


def _sub(a, b):
    return (sp.expand(a[0] - b[0]), sp.expand(a[1] - b[1]))


def _add(a, b):
    return (sp.expand(a[0] + b[0]), sp.expand(a[1] + b[1]))


def _mul(t, a):
    return (sp.expand(t * a[0]), sp.expand(t * a[1]))


def _same(a, b):
    return _eq(a[0], b[0]) and _eq(a[1], b[1])


def _between(t):
    return _sgn(t) >= 0 and _sgn(1 - t) >= 0


def _segment_parameters(a, b, c, d):
    r = _sub(b, a)
    s = _sub(d, c)
    den = _cross(r, s)
    out = []
    if not _eq(den, 0):
        t = sp.cancel(_cross(_sub(c, a), s) / den)
        u = sp.cancel(_cross(_sub(c, a), r) / den)
        if _between(t) and _between(u):
            out.append((t, u))
        return out
    if not _eq(_cross(_sub(c, a), r), 0):
        return out
    rr = sp.expand(r[0] ** 2 + r[1] ** 2)
    if _eq(rr, 0):
        return out
    for p in (c, d):
        t = sp.cancel(((p[0] - a[0]) * r[0] + (p[1] - a[1]) * r[1]) / rr)
        if _between(t):
            out.append((t, None))
    ss = sp.expand(s[0] ** 2 + s[1] ** 2)
    if not _eq(ss, 0):
        for p in (a, b):
            u = sp.cancel(((p[0] - c[0]) * s[0] + (p[1] - c[1]) * s[1]) / ss)
            if _between(u):
                out.append((None, u))
    return out


def _winding(point, paths):
    _, y = point
    w = 0
    for path in paths:
        for a, b in zip(path, path[1:] + path[:1], strict=True):
            ay = _sgn(a[1] - y)
            by = _sgn(b[1] - y)
            orient = _sgn(_cross(_sub(b, a), _sub(point, a)))
            if ay <= 0 < by and orient > 0:
                w += 1
            elif by <= 0 < ay and orient < 0:
                w -= 1
    return w


def _filled(w, rule):
    return {
        "even_odd": w % 2 != 0,
        "nonzero": w != 0,
        "positive": w > 0,
        "negative": w < 0,
        "winding_at_least_two": abs(w) >= 2,
    }[rule]


def _sample_side(a, b, left, all_edges):
    mid = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
    r = _sub(b, a)
    n = (-r[1], r[0])
    if not left:
        n = (-n[0], -n[1])
    # Find the first positive intersection of the normal ray with another
    # atomic segment and sample strictly before it.  This is exact and avoids
    # a scale-dependent epsilon.
    positive = []
    for c, d in all_edges:
        rr = n
        ss = _sub(d, c)
        den = _cross(rr, ss)
        if _eq(den, 0):
            continue
        t = sp.cancel(_cross(_sub(c, mid), ss) / den)
        u = sp.cancel(_cross(_sub(c, mid), rr) / den)
        if _sgn(t) > 0 and _between(u):
            positive.append(t)
    if positive:
        first = min(positive, key=cmp_to_key(lambda x, y: compare_exact_reals(x, y)))
        eps = sp.cancel(first / 2)
    else:
        eps = sp.Rational(1, 2)
    return _add(mid, _mul(eps, n))


def _area(cycle):
    return sp.expand(
        sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(cycle, cycle[1:] + cycle[:1], strict=True))
        / 2
    )


def _point_in_cycle(p, cycle):
    return _winding(p, (cycle,)) != 0


def polygonal_region_from_paths(paths, *, fill_rule="nonzero"):
    """Construct the exact planar filled set selected by a winding fill rule."""
    if fill_rule not in {"even_odd", "nonzero", "positive", "negative", "winding_at_least_two"}:
        raise ValueError("unsupported fill_rule")
    norm = []
    for raw in paths:
        p = [tuple(map(sp.sympify, q)) for q in raw]
        if len(p) > 1 and _same(p[0], p[-1]):
            p.pop()
        if len(p) < 2 or any(len(q) != 2 for q in p):
            raise ValueError("paths must contain planar points")
        norm.append(tuple(p))
    segments = [
        (a, b) for p in norm for a, b in zip(p, p[1:] + p[:1], strict=True) if not _same(a, b)
    ]
    params = [{sp.Integer(0), sp.Integer(1)} for _ in segments]
    for i, (a, b) in enumerate(segments):
        for j in range(i + 1, len(segments)):
            c, d = segments[j]
            for t, u in _segment_parameters(a, b, c, d):
                if t is not None:
                    params[i].add(sp.cancel(t))
                if u is not None:
                    params[j].add(sp.cancel(u))
    atoms = []
    for (a, b), ts in zip(segments, params, strict=True):
        r = _sub(b, a)
        ordered = sorted(ts, key=cmp_to_key(lambda x, y: compare_exact_reals(x, y)))
        for s, t in zip(ordered, ordered[1:], strict=False):
            if not _eq(s, t):
                atoms.append((_add(a, _mul(s, r)), _add(a, _mul(t, r))))
    # Deduplicate geometric atomic segments; winding is evaluated against original paths so overlap multiplicity is retained.
    unique = []
    for a, b in atoms:
        if not any(
            (_same(a, c) and _same(b, d)) or (_same(a, d) and _same(b, c)) for c, d in unique
        ):
            unique.append((a, b))
    boundary = []
    for a, b in unique:
        lf = _filled(_winding(_sample_side(a, b, True, unique), norm), fill_rule)
        rf = _filled(_winding(_sample_side(a, b, False, unique), norm), fill_rule)
        if lf != rf:
            boundary.append((a, b) if lf else (b, a))
    outgoing = defaultdict(list)

    def key(p):
        return tuple(sp.srepr(x) for x in p)

    points = {key(p): p for e in boundary for p in e}
    for a, b in boundary:
        outgoing[key(a)].append(key(b))
    cycles = []
    while outgoing:
        start = min(outgoing)
        cur = start
        cyc = []
        while True:
            cyc.append(points[cur])
            choices = outgoing[cur]
            nxt = min(choices)
            choices.remove(nxt)
            if not choices:
                del outgoing[cur]
            cur = nxt
            if cur == start:
                break
        if len(cyc) >= 3:
            cycles.append(tuple(cyc))
    outers = [c for c in cycles if _sgn(_area(c)) > 0]
    holes = [c for c in cycles if _sgn(_area(c)) < 0]
    components = []
    for outer in outers:
        owned = []
        for hole in holes:
            if _point_in_cycle(hole[0], outer):
                # choose smallest containing outer by exact absolute area
                containers = [o for o in outers if _point_in_cycle(hole[0], o)]
                chosen = min(
                    containers,
                    key=cmp_to_key(lambda a, b: compare_exact_reals(abs(_area(a)), abs(_area(b)))),
                )
                if chosen is outer:
                    owned.append(Polygon(hole))
        components.append(PolygonalComponent(Polygon(outer), owned))
    return PolygonalSet(components)
