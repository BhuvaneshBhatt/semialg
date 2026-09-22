"""Contracts for the generated function-by-function root API test adequacy audit."""

from __future__ import annotations

import inspect
import tomllib
from pathlib import Path

import semialg
from scripts.generate_root_api_test_adequacy import render

ROOT = Path(__file__).resolve().parents[1]


def test_test_adequacy_registry_covers_exactly_root_functions() -> None:
    data = tomllib.loads((ROOT / "tests/root_api_test_adequacy.toml").read_text())["apis"]
    public = {n for n in semialg.__all__ if inspect.isfunction(getattr(semialg, n))}
    assert set(data) == public


def test_test_adequacy_registry_is_reproducible() -> None:
    toml_text, report = render()
    assert (ROOT / "tests/root_api_test_adequacy.toml").read_text() == toml_text
    assert (ROOT / "docs/quality/root-api-test-adequacy.md").read_text() == report


def test_test_adequacy_status_matches_missing_dimensions() -> None:
    data = tomllib.loads((ROOT / "tests/root_api_test_adequacy.toml").read_text())["apis"]
    for item in data.values():
        assert (item["status"] == "adequate") == (item["missing"] == [])
        assert "nominal" in item["dimensions"]
