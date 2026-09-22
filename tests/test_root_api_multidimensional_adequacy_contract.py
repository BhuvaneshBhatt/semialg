"""Enforce semantic dimensions, not merely call counts, for deepened root APIs."""

from __future__ import annotations

import ast
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "tests/root_api_adequacy_dimensions.toml"
ALLOWED = {
    "nominal",
    "metamorphic",
    "boundary-degenerate",
    "independent-oracle",
    "independent-semantic",
    "round-trip-presentation",
}


def _registry():
    return tomllib.loads(REGISTRY.read_text())["apis"]


def _test_directly_calls(ref: str, api: str) -> bool:
    file_name, test_name = ref.split("::", 1)
    tree = ast.parse((ROOT / file_name).read_text(), filename=file_name)
    fn = next(
        n
        for n in tree.body
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name == test_name
    )
    for node in ast.walk(fn):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Name) and node.func.id == api:
            return True
        if (
            isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "sa"
            and node.func.attr == api
        ):
            return True
    return False


def test_multidimensional_registry_covers_all_previously_minimum_depth_apis():
    registry = _registry()
    assert len(registry) == 79
    for api, item in registry.items():
        dimensions = set(item["dimensions"])
        assert "nominal" in dimensions
        assert len(dimensions) >= 2, (api, dimensions)
        assert dimensions <= ALLOWED
        assert item["evidence"]
        assert all(_test_directly_calls(ref, api) for ref in item["evidence"])


def test_multidimensional_registry_is_reproducible():
    import subprocess
    import sys

    before = REGISTRY.read_text()
    subprocess.run(
        [sys.executable, str(ROOT / "scripts/generate_root_api_adequacy_dimensions.py")],
        check=True,
        cwd=ROOT,
    )
    assert REGISTRY.read_text() == before
