from __future__ import annotations

from pathlib import Path

import pytest

SEMANTIC_TIERS = ("contract", "documentation", "integration", "property", "wheel")
_DOCUMENTATION_FILES = {
    "test_demo_notebook.py",
    "test_documentation_code_blocks.py",
    "test_documentation_example_gallery.py",
    "test_documentation_information_architecture.py",
    "test_documentation_semantic_examples.py",
    "test_documented_defaults_contract.py",
    "test_examples_executable.py",
    "test_mkdocs_navigation_coverage.py",
    "test_public_api_documentation.py",
    "test_region_documentation_contracts.py",
    "test_semialgebraic_geometry_documentation.py",
}


def _semantic_tier(path: Path) -> str:
    parts = path.parts
    if "release" in parts:
        return "wheel"
    if "properties" in parts:
        return "property"
    if path.name in _DOCUMENTATION_FILES or "documentation" in path.name:
        return "documentation"
    if "contract" in path.name:
        return "contract"
    return "integration"


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    """Give every test one stable semantic tier independent of its filename shard."""

    for item in items:
        tiers = [tier for tier in SEMANTIC_TIERS if tier in item.keywords]
        if not tiers:
            tier = _semantic_tier(Path(str(item.path)))
            item.add_marker(getattr(pytest.mark, tier))
            tiers = [tier]
        if len(tiers) != 1:
            raise pytest.UsageError(f"{item.nodeid} has conflicting semantic tiers: {tiers}")
