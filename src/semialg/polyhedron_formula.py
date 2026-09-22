"""Exact implicit semantics for shell-based polyhedra."""

from __future__ import annotations

import sympy as sp

from .boundary_topology import PolyhedralShell, Polyhedron
from .derived_geometry import connected_components, is_bounded


def _face_triangles(shell: PolyhedralShell):
    for face in shell.faces:
        a = face[0]
        for i in range(1, len(face) - 1):
            yield (shell.vertices[a], shell.vertices[face[i]], shell.vertices[face[i + 1]])


def _triangle_formula(triangle, variables):
    a, b, c = (sp.Matrix(p) for p in triangle)
    x = sp.Matrix(variables)
    n = (b - a).cross(c - a)
    plane = sp.Eq(sp.expand(n.dot(x - a)), 0)
    tests = []
    for u, v in ((a, b), (b, c), (c, a)):
        tests.append(sp.Ge(sp.expand(n.dot((v - u).cross(x - u))), 0))
    return sp.And(plane, *tests)


def _shell_boundary_formula(shell, variables):
    return sp.Or(*(_triangle_formula(t, variables) for t in _face_triangles(shell)))


def _shell_solid_formula(shell, variables):
    boundary = _shell_boundary_formula(shell, variables)
    try:
        from .canonicalization import _shell_is_convex_boundary

        if _shell_is_convex_boundary(shell):
            from .standard_regions import Polytope

            return Polytope(shell.vertices).h_representation().as_formula(variables), boundary
    except (ValueError, TypeError, NotImplementedError):
        pass
    complement = sp.Not(boundary)
    bounded = tuple(
        c for c in connected_components(complement, variables) if is_bounded(c, variables)
    )
    return sp.Or(boundary, *bounded), boundary


def polyhedron_formula(polyhedron: Polyhedron, variables):
    """Return an exact formula for a shell-based polyhedron, including cavities.

    Each closed shell is interpreted by the bounded connected components of
    the complement of its exact polygonal boundary.  This avoids convexifying
    nonconvex shells or inventing a tetrahedral fill.
    """
    if not isinstance(polyhedron, Polyhedron):
        raise TypeError("polyhedron must be a boundary_topology.Polyhedron")
    variables = tuple(variables)
    if len(variables) != 3:
        raise ValueError("Polyhedron formulas require three variables")
    solids = []
    for component in polyhedron.components:
        outer, _ = _shell_solid_formula(component.outer_shell, variables)
        body = outer
        for cavity in component.cavity_shells:
            cavity_solid, cavity_boundary = _shell_solid_formula(cavity, variables)
            body = sp.And(body, sp.Or(sp.Not(cavity_solid), cavity_boundary))
        solids.append(body)
    return sp.Or(*solids) if solids else sp.false
