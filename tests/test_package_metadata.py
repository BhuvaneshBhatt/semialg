from __future__ import annotations

from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import semialg


def test_package_version_matches_project_metadata():
    root = Path(__file__).resolve().parents[1]
    with (root / "pyproject.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    assert semialg.__version__ == metadata["project"]["version"]
