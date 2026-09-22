"""Structural contracts for public imports and exception boundaries."""

from __future__ import annotations

import ast
from pathlib import Path

import semialg
from semialg import cad_region, region_coercion, region_predicates, symbolic_regions
from semialg._public_api import PUBLIC_EXPORTS
from semialg.cad_algorithms import meshing, numerical_boundaries

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "semialg"


def test_cad_region_exports_implementation_objects():
    assert cad_region.CADMesh is meshing.CADMesh
    assert cad_region.triangulate_cad_cell is meshing.triangulate_cad_cell
    assert cad_region.triangulate_cad_cells is meshing.triangulate_cad_cells
    assert cad_region.triangulate_cad_region is meshing.triangulate_cad_region
    assert cad_region.evaluate_delineable_curve is numerical_boundaries.evaluate_delineable_curve
    assert (
        cad_region.evaluate_delineable_surfaces is numerical_boundaries.evaluate_delineable_surfaces
    )


def test_symbolic_region_exports_coercion_and_predicates():
    assert symbolic_regions.as_semialgebraic_region is region_coercion.as_semialgebraic_region
    assert symbolic_regions.RegionElement is region_predicates.RegionElement
    assert symbolic_regions.RegionNotElement is region_predicates.RegionNotElement
    assert symbolic_regions.RegionSubset is region_predicates.RegionSubset
    assert symbolic_regions.RegionDisjoint is region_predicates.RegionDisjoint
    assert symbolic_regions.RegionEqual is region_predicates.RegionEqual


def test_top_level_public_namespace_has_one_registry_source_of_truth():
    assert semialg.__all__ == ["__version__", *PUBLIC_EXPORTS]
    assert len(PUBLIC_EXPORTS) == len(set(PUBLIC_EXPORTS))
    for name, module_name in PUBLIC_EXPORTS.items():
        assert module_name.startswith(".")
        assert getattr(semialg, name) is not None


def _is_broad_exception(handler: ast.ExceptHandler) -> bool:
    node = handler.type
    return isinstance(node, ast.Name) and node.id in {"Exception", "BaseException"}


def test_broad_exception_catches_are_confined_to_documented_boundaries():
    # Broad catches are limited to boundaries that execute arbitrary external
    # or user-defined behavior.
    allowed = {
        "conditional.py",
        "solution_geometry.py",
        "solve/integer/engine.py",
        "validation/checkers.py",
        "validation/runner.py",
    }
    offenders: list[tuple[str, int]] = []
    for path in SRC.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        rel = path.relative_to(SRC).as_posix()
        for node in ast.walk(tree):
            if isinstance(node, ast.ExceptHandler) and _is_broad_exception(node):
                if rel not in allowed:
                    offenders.append((rel, node.lineno))
    assert offenders == []
