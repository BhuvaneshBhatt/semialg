"""Deterministic exact canonicalization of polygonal/polyhedral geometry."""

from __future__ import annotations

import itertools

import sympy as sp

from ._zero_testing import certified_equal, certified_sign
from .boundary_topology import PolygonalSet, PolyhedralComponent, PolyhedralShell, Polyhedron
from .standard_regions import (
    Box,
    Hexahedron,
    Parallelepiped,
    Parallelogram,
    Polygon,
    Polytope,
    Simplex,
    TetrahedralComplex,
    Zonotope,
)


def _eq(a, b):
    return certified_equal(sp.sympify(a), sp.sympify(b)) is True


def _same(a, b):
    return len(a) == len(b) and all(_eq(x, y) for x, y in zip(a, b, strict=True))


def _cross(a, b, c):
    return sp.expand((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def _norm(vertices):
    pts = [tuple(map(sp.sympify, p)) for p in vertices]
    if len(pts) > 1 and _same(pts[0], pts[-1]):
        pts.pop()
    changed = True
    while changed and len(pts) > 3:
        changed = False
        for i in range(len(pts)):
            if _eq(_cross(pts[i - 1], pts[i], pts[(i + 1) % len(pts)]), 0):
                pts.pop(i)
                changed = True
                break
    if len(pts) < 3:
        raise ValueError("polygon must have at least three non-collinear vertices")
    area = sp.expand(
        sum(
            pts[i][0] * pts[(i + 1) % len(pts)][1] - pts[(i + 1) % len(pts)][0] * pts[i][1]
            for i in range(len(pts))
        )
    )
    if certified_sign(area) == -1:
        pts.reverse()
    start = min(range(len(pts)), key=lambda i: tuple(sp.default_sort_key(x) for x in pts[i]))
    return tuple(pts[start:] + pts[:start])


def _axis_box(pts):
    d = len(pts[0])
    coords = [tuple(sorted({p[j] for p in pts}, key=sp.default_sort_key)) for j in range(d)]
    if all(len(c) == 2 for c in coords) and len(pts) == 2**d:
        expected = {tuple(v) for v in itertools.product(*coords)}
        if set(pts) == expected:
            return Box(tuple((c[0], c[1]) for c in coords))
    return None


def canonicalize_polygon(region):
    """Return the strongest certified canonical representation of polygonal geometry."""
    if isinstance(region, PolygonalSet):
        return region
    pts = _norm(region.vertices if isinstance(region, Polygon) else region)
    if len(pts) == 3:
        return Simplex(pts)
    box = _axis_box(pts)
    if box is not None:
        return box
    if len(pts) == 4 and all(_eq(pts[0][j] + pts[2][j], pts[1][j] + pts[3][j]) for j in range(2)):
        return Parallelogram(
            pts[0],
            (
                tuple(pts[1][j] - pts[0][j] for j in range(2)),
                tuple(pts[3][j] - pts[0][j] for j in range(2)),
            ),
        )
    return Polygon(pts)


def _extreme_polytope(p):
    if p.dimension() != p.ambient_dimension():
        return p
    try:
        incidence = p.incidence()
        extreme = tuple(
            v for v, facets in zip(p.vertices, incidence.vertex_facets, strict=True) if facets
        )
        return Polytope(extreme) if len(extreme) < len(p.vertices) else p
    except (ValueError, TypeError, NotImplementedError):
        return p


def _parallelepiped(p):
    d = p.dimension()
    pts = p.vertices
    if d != p.ambient_dimension() or len(pts) != 2**d:
        return None
    box = _axis_box(pts)
    if box is not None:
        return box
    lattice = p.face_lattice()
    edges = lattice.faces_of_dimension(1)
    adjacency = {i: set() for i in range(len(pts))}
    for e in edges:
        a, b = e.vertex_indices
        adjacency[a].add(b)
        adjacency[b].add(a)
    origin = min(range(len(pts)), key=lambda i: tuple(sp.default_sort_key(x) for x in pts[i]))
    neighbors = sorted(adjacency[origin])
    if len(neighbors) != d:
        return None
    vecs = [tuple(sp.expand(pts[n][j] - pts[origin][j]) for j in range(d)) for n in neighbors]
    expected = []
    for bits in itertools.product((0, 1), repeat=d):
        expected.append(
            tuple(
                sp.expand(pts[origin][j] + sum(bits[k] * vecs[k][j] for k in range(d)))
                for j in range(d)
            )
        )
    if all(any(_same(q, p0) for p0 in pts) for q in expected):
        return Parallelepiped(pts[origin], vecs)
    return None


def _zonotope(p):
    """Recognize a vertex-presented zonotope when an exact generator replay succeeds."""
    d = p.dimension()
    if d <= 0 or d != p.ambient_dimension():
        return None
    lattice = p.face_lattice()
    edges = lattice.faces_of_dimension(1)
    from .exact_arithmetic import compare_exact_reals

    generators = []
    for edge in edges:
        i, j = edge.vertex_indices
        vec = tuple(
            sp.expand(p.vertices[j][k] - p.vertices[i][k]) for k in range(p.ambient_dimension())
        )
        first = next((x for x in vec if not _eq(x, 0)), None)
        if first is None:
            continue
        if compare_exact_reals(first, 0) < 0:
            vec = tuple(-x for x in vec)
        if not any(_same(vec, g) for g in generators):
            generators.append(vec)
    if not generators or len(generators) > 10:
        return None
    center = tuple(
        sp.cancel(sum(v[j] for v in p.vertices) / len(p.vertices))
        for j in range(p.ambient_dimension())
    )
    origin = tuple(
        sp.cancel(center[j] - sum(g[j] for g in generators) / 2)
        for j in range(p.ambient_dimension())
    )
    generated = [
        tuple(
            sp.expand(origin[j] + sum(bits[k] * generators[k][j] for k in range(len(generators))))
            for j in range(p.ambient_dimension())
        )
        for bits in itertools.product((0, 1), repeat=len(generators))
    ]
    try:
        replay = _extreme_polytope(Polytope(generated))
    except (ValueError, TypeError, NotImplementedError):
        return None
    if len(replay.vertices) == len(p.vertices) and all(
        any(_same(v, w) for w in p.vertices) for v in replay.vertices
    ):
        return Zonotope(origin, generators)
    return None


def _canonical_polytope(p):
    p = _extreme_polytope(p)
    d = p.dimension()
    n = len(p.vertices)
    if n == d + 1:
        return Simplex(p.vertices)
    para = _parallelepiped(p)
    if para is not None:
        return para
    zono = _zonotope(p)
    if zono is not None:
        return zono
    if d == 3 and n == 8:
        try:
            return Hexahedron(p.vertices)
        except ValueError:
            pass
    return p


def _merge_coplanar_faces(shell):
    """Merge adjacent coplanar polygonal faces without changing shell geometry."""
    faces = [tuple(f) for f in shell.faces]
    verts = shell.vertices

    def plane(face):
        a = sp.Matrix(verts[face[0]])
        for i in range(1, len(face) - 1):
            n = (sp.Matrix(verts[face[i]]) - a).cross(sp.Matrix(verts[face[i + 1]]) - a)
            if any(not _eq(v, 0) for v in n):
                return a, n
        return None, None

    def coplanar(f, g):
        a, n = plane(f)
        if n is None:
            return False
        return all(_eq(n.dot(sp.Matrix(verts[i]) - a), 0) for i in g)

    def merge(f, g):
        for i in range(len(f)):
            a, b = f[i], f[(i + 1) % len(f)]
            for j in range(len(g)):
                if g[j] == b and g[(j + 1) % len(g)] == a:
                    pf = []
                    k = (i + 1) % len(f)
                    while True:
                        pf.append(f[k])
                        k = (k + 1) % len(f)
                        if k == (i + 1) % len(f):
                            break
                    pg = []
                    k = (j + 1) % len(g)
                    while True:
                        pg.append(g[k])
                        k = (k + 1) % len(g)
                        if k == (j + 1) % len(g):
                            break
                    # Paths begin at b/a respectively; omit the shared-edge
                    # endpoints duplicated by concatenation.
                    cycle = pf[:-1] + pg[:-1]
                    cleaned = []
                    for v in cycle:
                        if not cleaned or v != cleaned[-1]:
                            cleaned.append(v)
                    return tuple(cleaned)
        return None

    changed = True
    while changed:
        changed = False
        for i in range(len(faces)):
            for j in range(i + 1, len(faces)):
                if coplanar(faces[i], faces[j]):
                    joined = merge(faces[i], faces[j])
                    if joined is not None and len(set(joined)) >= 3:
                        faces = [f for k, f in enumerate(faces) if k not in (i, j)] + [joined]
                        changed = True
                        break
            if changed:
                break
    try:
        return PolyhedralShell(verts, faces, deduplicate_vertices=False)
    except (ValueError, TypeError):
        return shell


def _shell_is_convex_boundary(shell):
    try:
        hull = Polytope(shell.vertices)
        return {frozenset(f.vertex_indices) for f in hull.facets()} == {
            frozenset(f) for f in shell.faces
        }
    except (ValueError, TypeError, NotImplementedError):
        return False


def canonicalize_polyhedron(region):
    """Return a certified canonical polyhedral representation without unsafe convexification."""
    if isinstance(region, Polytope):
        return _canonical_polytope(region)
    if isinstance(region, TetrahedralComplex):
        if len(region.tetrahedra) == 1:
            return Simplex(region.tetrahedra[0].vertices)
        return region
    if isinstance(region, PolyhedralShell):
        return (
            _canonical_polytope(Polytope(region.vertices))
            if _shell_is_convex_boundary(region)
            else region
        )
    if isinstance(region, Polyhedron):
        if len(region.components) == 1 and not region.components[0].cavity_shells:
            shell = region.components[0].outer_shell
            if _shell_is_convex_boundary(shell):
                return _canonical_polytope(Polytope(shell.vertices))
        # Reclassify all supplied shells by geometric nesting.  This repairs
        # callers that supplied nested closed shells as separate components.
        shells = []
        for component in region.components:
            shells.append(component.outer_shell)
            shells.extend(component.cavity_shells)
        if len(shells) > 1:
            try:
                from .derived_geometry import is_subset
                from .polyhedron_formula import _shell_boundary_formula, _shell_solid_formula

                vars_ = sp.symbols("x0:3", real=True)
                solids = [_shell_solid_formula(sh, vars_)[0] for sh in shells]
                boundaries = [_shell_boundary_formula(sh, vars_) for sh in shells]
                contains = [[False] * len(shells) for _ in shells]
                for i in range(len(shells)):
                    for j in range(len(shells)):
                        if i == j:
                            continue
                        if _shell_is_convex_boundary(shells[i]):
                            h = Polytope(shells[i].vertices).h_representation()
                            ok = True
                            for point in shells[j].vertices:
                                for row in range(h.matrix.rows):
                                    value = sp.expand(
                                        sum(h.matrix[row, k] * point[k] for k in range(3))
                                        - h.offsets[row, 0]
                                    )
                                    from .exact_arithmetic import compare_exact_reals

                                    if compare_exact_reals(value, 0) > 0:
                                        ok = False
                                        break
                                if not ok:
                                    break
                            contains[i][j] = ok
                        else:
                            contains[i][j] = is_subset(boundaries[j], solids[i], vars_)
                depth = [
                    sum(contains[i][j] for i in range(len(shells)) if i != j)
                    for j in range(len(shells))
                ]
                components = []
                for i, sh in enumerate(shells):
                    if depth[i] % 2:
                        continue
                    children = [
                        shells[j]
                        for j in range(len(shells))
                        if depth[j] == depth[i] + 1 and contains[i][j]
                    ]
                    components.append(PolyhedralComponent(sh, children))
                return Polyhedron(tuple(components))
            except (ValueError, TypeError, NotImplementedError):
                pass
        return region
    raise TypeError("expected polyhedral geometry")


def canonicalize_region(region):
    """Canonicalize a supported region using exact structural recognition."""
    if isinstance(region, (Polygon, PolygonalSet)):
        return canonicalize_polygon(region)
    if isinstance(region, (Polytope, TetrahedralComplex, PolyhedralShell, Polyhedron)):
        return canonicalize_polyhedron(region)
    return region
