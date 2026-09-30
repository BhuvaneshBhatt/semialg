"""Source distributions must retain reproducible documentation-quality inputs."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_documentation_quality_registries_and_assets_are_source_files():
    required = (
        "docs/reference/primary_api_manifest.toml",
        "docs/reference/root_function_documentation_adequacy.toml",
        "docs/assets/decision-flow.svg",
        "docs/assets/cad-flow.svg",
        "docs/assets/analysis-flow.svg",
        "scripts/generate_primary_api_manifest.py",
        "scripts/generate_root_api_documentation_adequacy.py",
    )
    assert [name for name in required if not (ROOT / name).is_file()] == []


def test_manifest_includes_documentation_registries_and_assets():
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "recursive-include docs *.md *.toml *.svg" in manifest
