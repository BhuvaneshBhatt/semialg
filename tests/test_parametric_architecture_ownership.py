"""Permanent ownership contracts for the consolidated parametric-CAD subsystem."""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "semialg"


def _top_level_relative_imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imports: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.level and node.module:
            imports.add(node.module)
    return imports


def test_obsolete_parametric_ownership_modules_are_absent() -> None:
    assert not (SRC / "generic").exists()
    assert not (SRC / "parameter_stratification.py").exists()


def test_parametric_support_has_one_way_ownership() -> None:
    parametric = SRC / "decomposition" / "parametric.py"
    support = SRC / "decomposition" / "_parametric_support.py"
    assert support.exists()
    imports = _top_level_relative_imports(support)
    # Support may depend on CAD primitives, formulas, and topology, but must not
    # import the high-level parametric façade or parameter-feasibility API.
    assert "parametric" not in imports
    assert "parameters" not in imports
    assert "parameter_stratification" not in imports
    assert "generic" not in imports
    assert "._parametric_support" in parametric.read_text(encoding="utf-8")


def test_region_coercion_has_single_implementation_owner() -> None:
    implementations = []
    for path in SRC.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == "as_semialgebraic_region":
                implementations.append(path.relative_to(SRC).as_posix())
    # region_predicates intentionally has a lazy import adapter to break the
    # symbolic-region import cycle; region_coercion owns the real conversion.
    assert sorted(implementations) == ["region_coercion.py", "region_predicates.py"]
    predicate_source = (SRC / "region_predicates.py").read_text(encoding="utf-8")
    assert (
        "from .region_coercion import as_semialgebraic_region as coerce_region" in predicate_source
    )


def test_assignment_conversion_is_shared_not_reimplemented() -> None:
    zero_dim = (SRC / "solve" / "zero_dimensional.py").read_text(encoding="utf-8")
    rur = (SRC / "algebraic" / "rational_univariate" / "representation.py").read_text(
        encoding="utf-8"
    )
    assert "assignments_from_points" in zero_dim
    assert "assignments_from_points" in rur
