from __future__ import annotations

from pathlib import Path

import tomllib

import semialg


def test_package_version_matches_project_metadata():
    root = Path(__file__).resolve().parents[1]
    with (root / "pyproject.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    assert semialg.__version__ == metadata["project"]["version"]


def test_package_version_matches_changelog_release_heading():
    root = Path(__file__).resolve().parents[1]
    changelog = (root / "CHANGELOG.md").read_text(encoding="utf-8")
    assert f"## {semialg.__version__}" in changelog
