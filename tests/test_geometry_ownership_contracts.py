from __future__ import annotations

import ast
from pathlib import Path

import sympy as sp

import semialg
from semialg import RegionElement
from semialg.standard_regions import (
    AffineHalfSpace,
    AffineSpace,
    Ball,
    BooleanRegion,
    Box,
    Capsule,
    Cone,
    Cylinder,
    Ellipsoid,
    EllipsoidBoundary,
    FilledTorus,
    HalfSpace,
    Hyperplane,
    Interval,
    Parallelepiped,
    Parallelogram,
    ParametricRegion,
    Point,
    Polygon,
    PolyhedralCone,
    Polytope,
    Prism,
    Pyramid,
    Ray,
    Simplex,
    Sphere,
    SphericalShell,
    Stadium,
    TetrahedralComplex,
    Torus,
)

ROOT = Path(__file__).resolve().parents[1]


def _regions():
    t = sp.symbols("t", real=True)
    triangle = Simplex(((0, 0), (1, 0), (0, 1)))
    tetra = Simplex(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)))
    return (
        (Point((1, 2)), (1, 2)),
        (AffineSpace((0, 0), ((1, 0),)), (2, 0)),
        (Hyperplane((1, 0), (0, 0)), (0, 2)),
        (HalfSpace((1, 0), (0, 0)), (-1, 0)),
        (Ray((0, 0), (1, 0)), (2, 0)),
        (AffineHalfSpace((0, 0), ((1, 0),), (0, 1)), (2, 3)),
        (Interval(0, 2), (1,)),
        (Box(((0, 1), (0, 2))), (sp.Rational(1, 2), 1)),
        (Simplex(((0, 0), (1, 0), (0, 1))), (sp.Rational(1, 4), sp.Rational(1, 4))),
        (Polygon(((0, 0), (2, 0), (2, 1), (0, 1))), (1, sp.Rational(1, 2))),
        (Polytope(((0, 0), (1, 0), (0, 1))), (sp.Rational(1, 4), sp.Rational(1, 4))),
        (PolyhedralCone((0, 0), rays=((1, 0), (0, 1))), (1, 1)),
        (tetra, (sp.Rational(1, 4),) * 3),
        (TetrahedralComplex((tetra,)), (sp.Rational(1, 4),) * 3),
        (Parallelogram((0, 0), ((1, 0), (0, 2))), (sp.Rational(1, 2), 1)),
        (Parallelepiped((0, 0, 0), ((1, 0, 0), (0, 1, 0), (0, 0, 1))), (sp.Rational(1, 2),) * 3),
        (Prism(triangle, (1, 1)), (sp.Rational(1, 2), sp.Rational(1, 2))),
        (Pyramid(Simplex(((0, 0), (1, 0))), (0, 1)), (sp.Rational(1, 3), sp.Rational(1, 3))),
        (Ball((0, 0), 2), (1, 0)),
        (Sphere((0, 0), 2), (2, 0)),
        (Ellipsoid((0, 0), ((4, 0), (0, 9))), (1, 0)),
        (EllipsoidBoundary((0, 0), ((4, 0), (0, 9))), (2, 0)),
        (SphericalShell((0, 0), (1, 2)), (sp.Rational(3, 2), 0)),
        (Cylinder((0, 0, 0), (0, 0, 2), 1), (0, 0, 1)),
        (Cone((0, 0, 0), (0, 0, 2), 1), (0, 0, 1)),
        (Torus((0, 0, 0), 2, 1), (3, 0, 0)),
        (FilledTorus((0, 0, 0), 2, 1), (2, 0, 0)),
        (Stadium((0, 0), (2, 0), 1), (1, 0)),
        (Capsule((0, 0, 0), (0, 0, 2), 1), (0, 0, 1)),
        (ParametricRegion((t,), ((t, 0, 1),), (t, t**2)), (sp.Rational(1, 2), sp.Rational(1, 4))),
        (BooleanRegion("union", (Ball((0, 0), 1), Ball((2, 0), 1))), (0, 0)),
    )


def test_standard_regions_lower_to_symbolic_regions():
    for region, _ in _regions():
        symbolic = region.as_semialgebraic_region()
        assert symbolic.ambient_dimension == region.ambient_dimension()
        assert len(symbolic.variables) == region.ambient_dimension()
        assert isinstance(region.as_formula(), sp.logic.boolalg.Boolean)


def test_membership_forms_agree():
    region = Ball((0, 0), 2)
    inside = (1, 0)
    outside = (3, 0)
    symbolic = region.as_semialgebraic_region()
    assert symbolic.contains(inside) is True
    assert symbolic.contains(outside) is False
    assert RegionElement(inside, region).evaluate() is True
    assert RegionElement(outside, region).evaluate() is False
    assert region.contains(inside) is True
    assert region.contains(outside) is False


def test_membership_api_is_public():
    assert semialg.RegionElement is RegionElement
    assert callable(semialg.contains_point)


def test_geometry_owns_formula_and_membership_entrypoints():
    path = ROOT / "src" / "semialg" / "standard_regions.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    overrides = []
    for node in tree.body:
        if not isinstance(node, ast.ClassDef) or node.name == "Geometry":
            continue
        methods = {item.name for item in node.body if isinstance(item, ast.FunctionDef)}
        overlap = methods & {"as_formula", "as_semialgebraic_region", "contains"}
        if overlap:
            overrides.append((node.name, sorted(overlap)))
    assert overrides == []


def test_geometry_module_has_no_eager_cad_dependency():
    path = ROOT / "src" / "semialg" / "standard_regions.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imported = []
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
        elif isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
    assert not any("cad" in name or "decision" in name for name in imported)
